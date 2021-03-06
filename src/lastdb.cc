// Copyright 2008, 2009, 2010, 2011, 2013, 2014 Martin C. Frith

// Read fasta-format sequences; construct a suffix array of them; and
// write the results to files.

#include "lastdb.hh"
#include "utilities.hh"
#include "SubsetSuffixArraySort.hh"
#include "DatabaseThread.hh"

using namespace std;

Alphabet alph;
LastdbArguments args;
DatabaseThread **dbThreads;
std::ifstream in;
DatabaseVolume *vol;
unsigned currentVolumeNumber = 0;

// Set up an alphabet (e.g. DNA or protein), based on the user options
void makeAlphabet(Alphabet &alph, const LastdbArguments &args) {
    if (!args.userAlphabet.empty()) alph.fromString(args.userAlphabet);
    else if (args.isProtein) alph.fromString(alph.protein);
    else alph.fromString(alph.dna);
}

// Does the first sequence look like it isn't really DNA?
bool isDubiousDna(const Alphabet &alph, const MultiSequence &multi) {
    const uchar *seq = multi.seqReader() + multi.seqBeg(0);
    unsigned dnaCount = 0;

    for (indexT i = 0; i < 100; ++i) {  // look at the first 100 letters
        uchar c = alph.canonical[seq[i]];
        if (c == alph.size) return false;  // we hit the end of the sequence early
        if (c < alph.size || c == alph.encode[(uchar) 'N']) ++dnaCount;
    }

    return dnaCount < 90;  // more than 10% unexpected letters
}

static void addSeeds(SubsetSuffixArray indexes[], unsigned &numOfIndexes,
                     const std::vector<std::string> &seedStrings,
                     const LastdbArguments &args, const Alphabet &alph) {
    for (unsigned x = 0; x < seedStrings.size(); ++x) {
        if (numOfIndexes >= maxNumOfIndexes) ERR("too many seed patterns");
        CyclicSubsetSeed &seed = indexes[numOfIndexes++].getSeed();
        seed.fromString(seedStrings[x], args.isCaseSensitive, alph.encode);
    }
}

// Set up the seed pattern(s), and return how many of them there are
unsigned makeSubsetSeeds(SubsetSuffixArray indexes[],
                         const LastdbArguments &args,
                         const Alphabet &alph) {
    unsigned numOfIndexes = 0;
    const std::string &a = alph.letters;

    for (unsigned x = 0; x < args.subsetSeedFiles.size(); ++x) {
        const std::string &name = args.subsetSeedFiles[x];
        std::vector<std::string> s = CyclicSubsetSeed::fromName(name);
        addSeeds(indexes, numOfIndexes, s, args, alph);
    }

    for (unsigned x = 0; x < args.seedPatterns.size(); ++x) {
        const std::string &mask = args.seedPatterns[x];
        std::vector<std::string> s = CyclicSubsetSeed::fromMask(a, mask);
        addSeeds(indexes, numOfIndexes, s, args, alph);
    }

    if (numOfIndexes == 0) {
        if (alph.letters == alph.dna) {
            const char *mask = "1T1001100101";  // YASS
            std::vector<std::string> s = CyclicSubsetSeed::fromMask(a, mask);
            addSeeds(indexes, numOfIndexes, s, args, alph);
        }
        else {
            std::vector<std::string> s = CyclicSubsetSeed::fromMask(a, "1");
            addSeeds(indexes, numOfIndexes, s, args, alph);
        }
    }

    return numOfIndexes;
}

void writePrjFile(const std::string &fileName, const LastdbArguments &args,
                  const Alphabet &alph, countT sequenceCount,
                  const std::vector<countT> &letterCounts,
                  unsigned volumes, unsigned numOfIndexes) {
    countT letterTotal = std::accumulate(letterCounts.begin(),
                                         letterCounts.end(), countT(0));

    std::ofstream f(fileName.c_str());
    f << "version=" <<

    #include "version.hh"
    << '\n';
    f << "alphabet=" << alph << '\n';
    f << "numofsequences=" << sequenceCount << '\n';
    f << "numofletters=" << letterTotal << '\n';
    f << "letterfreqs=";
    for (unsigned i = 0; i < letterCounts.size(); ++i) {
        if (i > 0) f << ' ';
        f << letterCounts[i];
    }
    f << '\n';

    if (!args.isCountsOnly) {
        f << "maxunsortedinterval=" << args.minSeedLimit << '\n';
        f << "masklowercase=" << args.isCaseSensitive << '\n';
        if (args.inputFormat != sequenceFormat::fasta)
            f << "sequenceformat=" << args.inputFormat << '\n';
        if (volumes + 1 > 0) f << "volumes=" << volumes << '\n';
        else f << "numofindexes=" << numOfIndexes << '\n';
    }

    if (!f) ERR("can't write file: " + fileName);
}

// Make one database volume, from one batch of sequences
void makeVolume(SubsetSuffixArray indexes[], unsigned numOfIndexes,
                const MultiSequence &multi, const LastdbArguments &args,
                const Alphabet &alph, const std::vector<countT> &letterCounts,
                const std::string &baseName) {
    LOG("writing prj and multi ...");
    //SEM_WAIT(io);
    writePrjFile(baseName + ".prj", args, alph, multi.finishedSequences(),
                 letterCounts, -1, numOfIndexes);
    multi.toFiles(baseName);
    //SEM_POST(io);

    SuffixArraySorter *s = new SuffixArraySorter();
    for (unsigned x = 0; x < numOfIndexes; ++x) {
        LOG("sorting...");
        indexes[x].sortIndex(multi.seqReader(), args.minSeedLimit, s);

        LOG("bucketing...");
        indexes[x].makeBuckets(multi.seqReader(), args.bucketDepth);

        LOG("writing suffix arrays...");
        indexT textLength = multi.finishedSize();
        //SEM_WAIT(io);
        if (numOfIndexes > 1) indexes[x].toFiles(baseName + char('a' + x), false, textLength);
        else indexes[x].toFiles(baseName, true, textLength);
        //SEM_POST(io);

        indexes[x].clearPositions();
    }
    delete s;
    LOG("done!");
    //SEM_POST(io);
}

// The max number of sequence letters, such that the total volume size
// is likely to be less than volumeSize bytes.  (This is crude, it
// neglects memory for the sequence names, and the fact that
// lowercase-masked letters and DNA "N"s aren't indexed.)
indexT maxLettersPerVolume(const LastdbArguments &args,
                           unsigned numOfIndexes) {
    std::size_t bytesPerLetter = isFastq(args.inputFormat) ? 2 : 1;
    std::size_t maxIndexBytesPerPosition = sizeof(indexT) + 1;
    maxIndexBytesPerPosition *= numOfIndexes;
    std::size_t x = bytesPerLetter * args.indexStep + maxIndexBytesPerPosition;
    std::size_t y = args.volumeSize / x * args.indexStep;
    indexT z = y;
    if (z < y) z = indexT(-1);
    return z;
}

// Read the next sequence, adding it to the MultiSequence and the SuffixArray
std::istream &
appendFromFasta(MultiSequence &multi,
                SubsetSuffixArray indexes[], unsigned numOfIndexes,
                const LastdbArguments &args, const Alphabet &alph,
                std::istream &in) {
    //SEM_WAIT(io);
    indexT maxSeqLen = maxLettersPerVolume(args, numOfIndexes);
    if (multi.finishedSequences() == 0) maxSeqLen = indexT(-1);

    std::size_t oldUnfinishedSize = multi.unfinishedSize();
    indexT oldFinishedSize = multi.finishedSize();

    if (args.inputFormat == sequenceFormat::fasta)
        multi.appendFromFastaLASTDB(in, maxSeqLen, args.unlimited);
    else
        multi.appendFromFastq(in, maxSeqLen);

    if (!multi.isFinished() && multi.finishedSequences() == 0)
        ERR("encountered a sequence that's too long");

    // encode the newly-read sequence
    alph.tr(multi.seqWriter() + oldUnfinishedSize,
            multi.seqWriter() + multi.unfinishedSize());

    if (isPhred(args.inputFormat))  // assumes one quality code per letter:
        checkQualityCodes(multi.qualityReader() + oldUnfinishedSize,
                          multi.qualityReader() + multi.unfinishedSize(),
                          qualityOffset(args.inputFormat));

    if (in && multi.isFinished() && !args.isCountsOnly)
        for (unsigned x = 0; x < numOfIndexes; ++x)
            indexes[x].addPositions(multi.seqReader(), oldFinishedSize,
                                    multi.finishedSize(), args.indexStep);
    //SEM_POST(io);

    return in;
}

void readPrjFile(const std::string &prjname,
                 const LastdbArguments &args,
                 const Alphabet &alph,
                 countT &sequenceCount,
                 std::vector<countT> &letterTotals,
                 unsigned &volumeNumber) {
    std::ifstream f(prjname.c_str());
    if (!f) ERR("can't open file: " + prjname);
    unsigned version = 0;

    std::string line, word;
    while (getline(f, line)) {
        std::istringstream iss(line);
        getline(iss, word, '=');
        if (word == "version") iss >> version;
        if (word == "numofsequences") iss >> sequenceCount;
        if (word == "letterfreqs") {
            for (int i = 0; i < alph.size; i++) {
                iss >> letterTotals[i];
            }
        }
        if (word == "volumes") iss >> volumeNumber;
    }

    if (version < 294 && version > 0)
        ERR("the fastdb files are old: please re-run fastdb");
}


void incrementalFormatWithNovel(const LastdbArguments &args,
                                const Alphabet &alph,
                                MultiSequence &multi,
                                SubsetSuffixArray indexes[maxNumOfIndexes],
                                const unsigned numOfIndexes,
                                const std::string &inputName) {
    unsigned volumeNumber = 0;
    countT sequenceCount = 0;
    std::vector<countT> letterCounts(alph.size);
    std::vector<countT> letterTotals(alph.size);

    // Read in the old Prj file
    readPrjFile(args.lastdbName + ".prj", args, alph,
                sequenceCount, letterTotals, volumeNumber);

    // Create a new set of volumes for the novel sequences
    std::ifstream inFileStream;
    std::istream &in = openIn(inputName, inFileStream);
    LOG("reading " << inputName << "...");

    while (in) {
        try {
            appendFromFasta(multi, indexes, numOfIndexes, args, alph, in);
            if (!args.isProtein &&
                args.userAlphabet.empty() &&
                sequenceCount == 0 &&
                isDubiousDna(alph, multi)) {
                std::cerr << "fastdb: that's some funny-lookin DNA\n";
            }

            if (multi.isFinished()) {
                ++sequenceCount;
                indexT lastSeq = multi.finishedSequences() - 1;
                alph.count(multi.seqReader() + multi.seqBeg(lastSeq),
                           multi.seqReader() + multi.seqEnd(lastSeq),
                           &letterCounts[0]);
                // memory-saving, which seems to be important on 32-bit systems:
                //!!if( args.isCountsOnly ) multi.reinitForAppending();
                if (args.isCountsOnly) multi.initForAppending(1);
            } else {
                const std::string baseName = args.lastdbName + stringify(volumeNumber++);
                makeVolume(indexes, numOfIndexes,
                           multi, args, alph, letterCounts, baseName);
                for (unsigned c = 0; c < alph.size; ++c) letterTotals[c] += letterCounts[c];
                letterCounts.assign(alph.size, 0);
                //!!multi.reinitForAppending();
                multi.initForAppending(1);
            }
        } catch (const std::exception &ex) {
            std::cerr << ex.what() << std::endl;
            std::cerr << "Encountered a malformed sequence. Ignoring sequence and continuing" << std::endl;
            multi.removeLatest();
        }
    }

    if (multi.finishedSequences() > 0) {
        if (volumeNumber == 0) {
            makeVolume(indexes, numOfIndexes,
                       multi, args, alph, letterCounts, args.lastdbName);
            return;
        }
        std::string baseName = args.lastdbName + stringify(volumeNumber++);
        makeVolume(indexes, numOfIndexes,
                   multi, args, alph, letterCounts, baseName);
    }

    for (unsigned c = 0; c < alph.size; ++c) letterTotals[c] += letterCounts[c];

    // Write out the updated Prj file
    writePrjFile(args.lastdbName + ".prj", args, alph,
                 sequenceCount, letterTotals, volumeNumber, numOfIndexes);
}

/*
void incrementalFormatWithNovel(const LastdbArguments &args,
                                const Alphabet &alph,
                                MultiSequence &multi,
                                SubsetSuffixArray indexes[maxNumOfIndexes],
                                const unsigned numOfIndexes,
                                const std::string &inputName)
{
	unsigned volumeNumber = 0;
	countT sequenceCount = 0;
	std::vector<countT> letterCounts( alph.size );
	std::vector<countT> letterTotals( alph.size );

	// Read in the old Prj file
	readPrjFile( args.lastdbName + ".prj", args, alph,
	             sequenceCount, letterTotals, volumeNumber);

	// Create a new set of volumes for the novel sequences
	std::ifstream inFileStream;
	std::istream& in = openIn( inputName, inFileStream );
	LOG( "reading " << inputName << "..." );

	while(in){
		try {
			appendFromFasta( multi, indexes, numOfIndexes, args, alph, in );
			if( !args.isProtein &&
			    args.userAlphabet.empty() &&
			    sequenceCount == 0 &&
			    isDubiousDna( alph, multi ) ){
				std::cerr << "fastdb: that's some funny-lookin DNA\n";
			}

			if( multi.isFinished() ){
				++sequenceCount;
				indexT lastSeq = multi.finishedSequences() - 1;
				alph.count( multi.seqReader() + multi.seqBeg(lastSeq),
				            multi.seqReader() + multi.seqEnd(lastSeq),
				            &letterCounts[0] );
				// memory-saving, which seems to be important on 32-bit systems:
				if( args.isCountsOnly ) multi.reinitForAppending();
			} else {
				const std::string baseName = args.lastdbName + stringify(volumeNumber++);
				makeVolume( indexes, numOfIndexes,
				            multi, args, alph, letterCounts, baseName );
				for( unsigned c = 0; c < alph.size; ++c ) letterTotals[c] += letterCounts[c];
				letterCounts.assign( alph.size, 0 );
				multi.reinitForAppending();
			}
		} catch (const std::exception &ex){
			std::cerr << ex.what() << std::endl;
			std::cerr << "Encountered a malformed sequence. Ignoring sequence and continuing" << std::endl;
			multi.removeLatest();
		}
	}

	if( multi.finishedSequences() > 0 ){
		if( volumeNumber == 0 ){
			makeVolume( indexes, numOfIndexes,
			            multi, args, alph, letterCounts, args.lastdbName );
			return;
		}
		std::string baseName = args.lastdbName + stringify(volumeNumber++);
		makeVolume( indexes, numOfIndexes,
		            multi, args, alph, letterCounts, baseName );
	}

	for( unsigned c = 0; c < alph.size; ++c ) letterTotals[c] += letterCounts[c];

	// Write out the updated Prj file
	writePrjFile( args.lastdbName + ".prj", args, alph,
	              sequenceCount, letterTotals, volumeNumber, numOfIndexes );
}
 */

void generateDifference(const std::string &filename, const string &dbname) {

    /*
     for (int i=0; i<volumes; i++){
            read in the volume sequence names into a set
            read through the file
            if the sequence is not found we output it.
     }
     */
}

void writeOuterPrj(unsigned numOfIndexes) {
    countT finalSequenceCount = 0;
    std::vector<countT> FinalLetterTotals(alph.size);
    for (int i = 0; i < args.threadNum; i++) {
        finalSequenceCount += dbThreads[i]->sequenceTotals;
        for (int j = 0; j < alph.size; j++) {
            FinalLetterTotals[j] += dbThreads[i]->letterTotals[j];
        }
    }

    countT letterTotal = std::accumulate(FinalLetterTotals.begin(),
                                         FinalLetterTotals.end(), countT(0));
    std::string fileName = args.lastdbName + ".prj";

    std::ofstream f(fileName.c_str());
    f << "version=" <<

    #include "version.hh"
    << '\n';
    f << "alphabet=" << alph << '\n';
    f << "numofsequences=" << finalSequenceCount << '\n';
    f << "numofletters=" << letterTotal << '\n';
    f << "letterfreqs=";
    for (unsigned i = 0; i < FinalLetterTotals.size(); ++i) {
        if (i > 0) f << ' ';
        f << FinalLetterTotals[i];
    }
    f << '\n';

    if (!args.isCountsOnly) {
        f << "maxunsortedinterval=" << args.minSeedLimit << '\n';
        f << "masklowercase=" << args.isCaseSensitive << '\n';
        if (args.inputFormat != sequenceFormat::fasta)
            f << "sequenceformat=" << args.inputFormat << '\n';
        if (currentVolumeNumber + 1 > 0) f << "volumes=" << (currentVolumeNumber + 1) << '\n';
        else f << "numofindexes=" << numOfIndexes << '\n';
    }

    if (!f) ERR("can't write file: " + fileName);
}

void createThreads() {
    dbThreads = new DatabaseThread *[args.threadNum];
    for (int i = 0; i < args.threadNum; i++) {
        dbThreads[i] = new DatabaseThread(alph.size, i);
    }
}

void lastdb(int argc, char **argv) {
    args.fromArgs(argc, argv);
    initializeSemaphores();
    makeAlphabet(alph, args);

    createThreads();

    alph.tr(dbThreads[0]->multi.seqWriter(),
            dbThreads[0]->multi.seqWriter() + dbThreads[0]->multi.unfinishedSize());

    unsigned numOfIndexes = dbThreads[0]->numOfIndexes;

    char defaultInputName[] = "-";
    char *defaultInput[] = {defaultInputName, 0};
    char **inputBegin = argv + args.inputStart;

    if (args.latestDatabase) {
        //!!
        // Run a diff across both databases, the formatted last and the new unformatted database
        // Write all of the sequences that dont appear in the formatted database to a temp file
        std::string tmpFile = generate_directory_name("");

        for (char **i = *inputBegin ? inputBegin : defaultInput; *i; ++i) {
            //generateDifference(tmpFile);
        }

        // Run the incrementalFormat with the new tmp file
        incrementalFormatWithNovel(args, alph, dbThreads[0]->multi,
                                   dbThreads[0]->indexes, numOfIndexes, tmpFile);

        // delete the old tmp file
        remove(tmpFile.c_str());
    } else if (args.novelSequenceFile) {
        for (char **i = *inputBegin ? inputBegin : defaultInput; *i; ++i)
            incrementalFormatWithNovel(args, alph, dbThreads[0]->multi,
                                       dbThreads[0]->indexes, numOfIndexes, string(*i));
    } else { // normal threading to build lastdb

        for (char **i = *inputBegin ? inputBegin : defaultInput; *i; ++i) {
            string inputName(*i);
            LOG("reading " << inputName << "...");
            in.open(inputName.c_str());

            //!! Create the DatabaseVolume for the first time here
            vol = new DatabaseVolume(args.inputFormat, args.lastdbName,
                                     currentVolumeNumber,
                                     numOfIndexes, alph.size);

            for (int i = 0; i < args.threadNum; i++) {
                LOG("Starting thread : " << i);
                dbThreads[i]->startThread();
            }

            for (int i = 0; i < args.threadNum; i++) {
                dbThreads[i]->joinThread();
                LOG("Joined thread : " << i);
            }
        }
    }

    /*
    //LOG("bucketing...");
    //!! Assuming we always produce only one index!
    std::cout << "indexTotal: " << vol->indexTotal << std::endl;
    //indexes[0].makeBuckets( multi.seqReader(), vol->indexTotal );
    indexes[0].buckets.v.resize( indexes[0].bucketSteps[0], vol->indexTotal );
    vol->writeBucketFile(indexes[0]);
    */

    // There may be a final volume which we never finished off.
    if (vol->indexTotal > 0) {
        // prj
        LOG("writing prj...");
        indexT textLength = vol->endsCoordinate;
        indexT bucketDepth = dbThreads[0]->indexes[0].defaultBucketDepth(vol->indexTotal);
        vol->prj->writePrjFile(args, dbThreads[0]->indexes[0],
                               textLength, vol->indexTotal,
                               bucketDepth);

        // buckets (bck)
        LOG("bucketing...");
        /*
        indexT bucketDepth = dbThreads[0]->indexes[0].defaultBucketDepth(vol->indexTotal);
        dbThreads[0]->indexes[0].makeBucketSteps( bucketDepth );
        dbThreads[0]->indexes[0].buckets.v.resize( dbThreads[0]->indexes[0].bucketSteps[0],
                                                   vol->indexTotal );
        LOG("writing buckets...");
        vol->writeBucketFile(dbThreads[0]->indexes[0]);
         */

        delete vol;
    }

    if (currentVolumeNumber > 1) {
        LOG("writing database prj...");
        writeOuterPrj(numOfIndexes);
    } else {
        renameFiles();
    }

    deleteThreads();
}

void renameFiles() {
    std::string n = args.lastdbName;
    rename((n + "0.prj").c_str(), (n + ".prj").c_str());
    rename((n + "0.ssp").c_str(), (n + ".ssp").c_str());
    rename((n + "0.tis").c_str(), (n + ".tis").c_str());
    rename((n + "0.sds").c_str(), (n + ".sds").c_str());
    rename((n + "0.des").c_str(), (n + ".des").c_str());
    rename((n + "0.suf").c_str(), (n + ".suf").c_str());
    rename((n + "0.bck").c_str(), (n + ".bck").c_str());
    if (args.inputFormat != sequenceFormat::fasta)
        rename((n + "0.qua").c_str(), (n + ".qua").c_str());
}

void deleteThreads() {
    for (int i = 0; i < args.threadNum; i++) {
        delete dbThreads[i];
    }
    delete[] dbThreads;
    destroySemaphores();
}

int main(int argc, char **argv) {
    try {
        lastdb(argc, argv);
        return EXIT_SUCCESS;
    } catch (const std::bad_alloc &e) {  // bad_alloc::what() may be unfriendly
        std::cerr << "fastdb: out of memory\n";
        return EXIT_FAILURE;
    } catch (const std::exception &e) {
        std::cerr << "fastdb: " << e.what() << '\n';
        return EXIT_FAILURE;
    } catch (int i) {
        return i;
    }
}
