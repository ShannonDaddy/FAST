// // Created by david on 16/04/16. //

#include "DatabaseThread.hh"
#include "SubsetSuffixArraySort.hh"

SEM_T io;

void DatabaseThread::accumulateAndFlushPrj(){
	LOG("accumulating prj...");
	vol->prj->accumulatePrj(letterCounts, sequenceCount);
	sequenceTotals += sequenceCount;
	//std::cout << "sequenceCount: " << sequenceCount << std::endl;
	//std::cout << "sequenceTotal: " << vol->prj->sequenceTotal << std::endl;
	for (unsigned c = 0; c < alph.size; ++c) letterTotals[c] += letterCounts[c];
	letterCounts.assign(alph.size, 0);
	sequenceCount = 0;
}

void DatabaseThread::accumulateAndFlushSuffixArrays(int x){
	LOG("accumulating suffix arrays...");
	indexT textLength = multi.finishedSize();
	//!! Original code snippet
	//if( numOfIndexes > 1 ) indexes[x].toFiles( baseName + char('a' + x), false, textLength );
	//else indexes[x].toFiles( baseName, true, textLength );

	//!! POUR OUT THE SUFFIX ARRAYS
	if (numOfIndexes > 1) vol->writePooledSubsetSuffixArray(indexes[x]);
	else vol->writePooledSubsetSuffixArray(indexes[x]);
	indexes[x].clearPositions();
}

void DatabaseThread::createSuffixArrays(int x){
	LOG("sorting...");
	indexes[x].sortIndex(multi.seqReader(), args.minSeedLimit, sorter);

	LOG("bucketing...");
	indexes[x].makeBuckets(multi.seqReader(), args.bucketDepth);
}

void DatabaseThread::accumulateAndFlushMulti(const LastdbArguments &args){
	LOG("accumulating multi...");
	vol->writePooledMultiSequence(multi, args);
	multi.reinitForAppending();
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
	// Check if the MultiSequence still has room in it for more
	if ( vol->isFinished() ) {
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

	} else {
		vol->prj->writePrjFile(args);
		replaceVolumeObject();
	}
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
				indexes[x].addPositions(multi.seqReader(), oldFinishedSize,
				                        multi.finishedSize(), args.indexStep);
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
		makeVolume(args, alph);
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
