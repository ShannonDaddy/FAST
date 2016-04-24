//
// Created by david on 16/04/16.
//

#include "DatabaseThread.hh"
#include "SubsetSuffixArraySort.hh"

SEM_T io;

void DatabaseThread::makeVolume( unsigned numOfIndexes,
                                 const LastdbArguments& args,
                                 const Alphabet& alph,
                                 const std::vector<countT>& letterCounts,
                                 const std::string& baseName )
{
	LOG( "writing prj and multi ..." );
	SEM_WAIT(io);
	writePrjFile( baseName + ".prj", args, alph, multi.finishedSequences(),
	              letterCounts, -1, numOfIndexes );
	multi.toFiles( baseName );
	SEM_POST(io);

	for( unsigned x = 0; x < numOfIndexes; ++x ){
		LOG( "sorting..." );
		indexes[x].sortIndex( multi.seqReader(), args.minSeedLimit, sorter );

		LOG( "bucketing..." );
		indexes[x].makeBuckets( multi.seqReader(), args.bucketDepth );

		LOG( "writing suffix arrays..." );
		indexT textLength = multi.finishedSize();
		SEM_WAIT(io);
		if( numOfIndexes > 1 ) indexes[x].toFiles( baseName + char('a' + x), false, textLength );
		else indexes[x].toFiles( baseName, true, textLength );
		SEM_POST(io);

		indexes[x].clearPositions();
	}
	LOG( "done!" );
	//SEM_POST(io);
}

std::istream&
DatabaseThread::readFasta( unsigned numOfIndexes,
                           const LastdbArguments& args,
                           const Alphabet& alph,
                           std::istream& in )
{
	indexT maxSeqLen = maxLettersPerVolume( args, numOfIndexes );
	if( multi.finishedSequences() == 0 ) maxSeqLen = indexT(-1);

	std::size_t seqCount = 0;
	while(seqCount < LOADSIZE) {
		SEM_WAIT(io);
		std::size_t oldUnfinishedSize = multi.unfinishedSize();
		indexT oldFinishedSize = multi.finishedSize();

		//SEM_WAIT(io);
		if (args.inputFormat == sequenceFormat::fasta)
			multi.appendFromFastaLASTDB(in, maxSeqLen, args.unlimited);
		else
			multi.appendFromFastq(in, maxSeqLen);
		seqCount++;
		//SEM_POST(io);

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

		if (!args.isProtein && args.userAlphabet.empty() &&
		    sequenceCount == 0 && isDubiousDna(alph, multi))
			std::cerr << "lastdb: that's some funny-lookin DNA\n";


		if (multi.isFinished()) {
			++sequenceCount;
			indexT lastSeq = multi.finishedSequences() - 1;
			alph.count(multi.seqReader() + multi.seqBeg(lastSeq),
			           multi.seqReader() + multi.seqEnd(lastSeq),
			           &letterCounts[0]);
			// memory-saving, which seems to be important on 32-bit systems:
			if (args.isCountsOnly) multi.reinitForAppending();
		} else prepareNextVolume();
		SEM_POST(io);
	}

	return in;
}

void DatabaseThread::prepareNextVolume()
{
	std::string baseName = args.lastdbName + stringify(volumeNumber++);
	makeVolume(numOfIndexes, args, alph,
	           letterCounts, baseName);
	for (unsigned c = 0; c < alph.size; ++c) letterTotals[c] += letterCounts[c];
	letterCounts.assign(alph.size, 0);
	multi.reinitForAppending();
}

void DatabaseThread::formatdb(const LastdbArguments &args,
                              const Alphabet &alph,
                              const unsigned numOfIndexes,
                              const std::string &inputName)
{
	LOG("reading " << inputName << "...");

	while (readFasta(numOfIndexes, args, alph, in))

		if( multi.finishedSequences() > 0 ){
			if( volumeNumber == 0 ) makeVolume(numOfIndexes, args, alph, letterCounts, args.lastdbName);
			else prepareNextVolume();
		}
}

void* DatabaseThread::threadEntry(void *args)
{
	((DatabaseThread *) args)->threadFunction();
	return NULL;
}

void DatabaseThread::threadFunction()
{
	formatdb(args, alph, numOfIndexes, currFile);
}

void DatabaseThread::startThread()
{
	pthread_create(&thread, NULL, threadEntry, this);
}

void DatabaseThread::joinThread()
{
	pthread_join(thread, NULL);
}

DatabaseThread::DatabaseThread(int s, int count):
//volumeNumber(0),
		volumeNumber(count),
		sequenceCount(0),
		rank(count)
{
	numOfIndexes = makeSubsetSeeds( indexes, args, alph );
	multi.initForAppending(1);
	letterCounts.resize(s);
	letterTotals.resize(s);
	sorter = new SuffixArraySorter();
}

DatabaseThread::~DatabaseThread()
{
	delete sorter;
}

void initializeSemaphores()
{
#ifdef __APPLE__
	sem_unlink("/io");
  if (( io = sem_open("/io", O_CREAT, 0644, 1)) == SEM_FAILED ) {
    perror("sem_open");
    exit(EXIT_FAILURE);
  }
#elif __linux
	sem_init(&io, 0, 1);
#endif
}

void destroySemaphores()
{
#ifdef __APPLE__
	sem_unlink("/io");
#elif __linux
	sem_destroy(&io);
#endif
}
