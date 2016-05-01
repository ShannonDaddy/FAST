// // Created by david on 16/04/16. //

#include "DatabaseThread.hh"
#include "SubsetSuffixArraySort.hh"

SEM_T io;

void DatabaseThread::accumulateAndFlushPrj(){
	LOG("accumulating prj...");
	vol->prj->accumulatePrj(letterCounts, sequenceCount);
	sequenceTotals += sequenceCount;

	for (unsigned c = 0; c < alph.size; ++c) letterTotals[c] += letterCounts[c];
	letterCounts.assign(alph.size, 0);
	sequenceCount = 0;
}

void DatabaseThread::accumulateAndFlushSuffixArrays(int x){
	LOG("accumulating suffix arrays...");
	//!! Original code snippet
	//if( numOfIndexes > 1 ) indexes[x].toFiles( baseName + char('a' + x), false, textLength );
	//else indexes[x].toFiles( baseName, true, textLength );

	vol->writePooledSubsetSuffixArray(indexes[x]);
	indexes[x].clearPositions();
}

void DatabaseThread::createSuffixArrays(int x){
	LOG("sorting...");
	indexes[x].sortIndex(multi.seqReader(), args.minSeedLimit, sorter);

	vol->indexTotal += indexes[x].index.size();
}

void DatabaseThread::accumulateAndFlushMulti(const LastdbArguments &args){
	LOG("accumulating multi...");

	unsigned endsLastCoordinate = multi.ends.back();
	unsigned nameEndsLastCoordinate = multi.nameEnds.back();

	// Update the ends with the global position
	for(int i=0; i<multi.ends.v.size(); i++){ // ssp file
		multi.ends.v[i] += vol->endsCoordinate;
	}

	// Update the nameEnds with the global position
	for(int i=0; i<multi.nameEnds.v.size(); i++){ // sds file
		multi.nameEnds.v[i] += vol->nameEndsCoordinate;
	}

	vol->writePooledMultiSequence(multi, args, endsLastCoordinate, nameEndsLastCoordinate);

	multi.initForAppending(1);
}

void DatabaseThread::replaceVolumeObject(){

	delete vol;	 // Get rid of the old one
	currentVolumeNumber++; // Increment the numbering scheme
	vol = new DatabaseVolume(args.inputFormat, args.lastdbName,
	                         currentVolumeNumber, numOfIndexes, alph.size);
}

void DatabaseThread::makeVolume( const LastdbArguments& args,
                                 const Alphabet& alph)
{
	accumulateAndFlushPrj();

	for (unsigned x = 0; x < numOfIndexes; ++x) {
		createSuffixArrays(x);
		SEM_WAIT(io);
		accumulateAndFlushSuffixArrays(x);
		SEM_POST(io);
	}

	SEM_WAIT(io);
	accumulateAndFlushMulti(args);
	SEM_POST(io);

	LOG("done with voluming batch");
}

std::istream&
DatabaseThread::readFasta( unsigned numOfIndexes,
                           const LastdbArguments& args,
                           const Alphabet& alph,
                           std::istream& in )
{
	indexT maxSeqLen = maxLettersPerVolume( args, numOfIndexes );
	if( multi.finishedSequences() == 0 ) maxSeqLen = indexT(-1);

	SEM_WAIT(io);
	while(sequenceCount < LOADSIZE && in) {
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

		if (in && multi.isFinished() && !args.isCountsOnly) {
			for (unsigned x = 0; x < numOfIndexes; ++x) {
				indexes[x].addPositions(multi.seqReader(),
				                        oldFinishedSize,
				                        multi.finishedSize(),
				                        args.indexStep);
			}
		}

		if( multi.finishedSequences() > 0 ) {
			++sequenceCount;
			indexT lastSeq = multi.finishedSequences() - 1;
			alph.count(multi.seqReader() + multi.seqBeg(lastSeq),
			           multi.seqReader() + multi.seqEnd(lastSeq),
			           &letterCounts[0]);
		}
	}
	SEM_POST(io);

	return in;
}

void DatabaseThread::formatdb(const LastdbArguments &args,
                              const Alphabet &alph)
{
	while (readFasta(numOfIndexes, args, alph, in)){
		// Check if the MultiSequence still has room in it for more
		if ( vol->isFinished(multi.ends.back(), multi.nameEnds.back()) ) {
			makeVolume(args, alph);
		} else {
			SEM_WAIT(io);
			//!! Assuming we always produce only one index!
			LOG("bucketing...");
			indexT bucketDepth = indexes[0].defaultBucketDepth(vol->indexTotal);
			indexes[0].makeBucketSteps( bucketDepth );
			indexes[0].buckets.v.resize( indexes[0].bucketSteps[0],
			                             vol->indexTotal );
			LOG("writing buckets...");
			vol->writeBucketFile(indexes[0]);

			//indexT textLength = multi.finishedSize();
			LOG("writing prj...");
			indexT textLength = vol->endsCoordinate;
			vol->prj->writePrjFile(args, indexes[0], textLength, vol->indexTotal);
			replaceVolumeObject();
			SEM_POST(io);

			makeVolume(args, alph);
		}
	}
}

void* DatabaseThread::threadEntry(void *_args)
{
	((DatabaseThread *) _args)->formatdb(args, alph);
	return NULL;
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
		sequenceCount(0),
		sequenceTotals(0),
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
