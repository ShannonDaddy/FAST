//
// Created by david on 16/04/16.
//

#ifndef THREADEDLAST_DATABASETHREAD_HH
#define THREADEDLAST_DATABASETHREAD_HH

#include "lastdb.hh"

class DatabaseThread{

		friend class DatabaseVolume;

//private:
public:
		pthread_t thread;

		MultiSequence multi;
		SubsetSuffixArray indexes[maxNumOfIndexes];

		unsigned numOfIndexes;
		// Constantly flushed for each volume chunks prj portion
		countT sequenceCount;
		// Keep track of all of the sequences used for the final prj file.
		countT sequenceTotals;
		// Constantly flushed for each volume chunks prj portion
		std::vector<countT> letterCounts;
		// Keep track of all of the letters used for the final prj file.
		std::vector<countT> letterTotals;

		int rank;
		SuffixArraySorter *sorter;

		void formatdb(const LastdbArguments &args,
		              const Alphabet &alph);

		void makeVolume(const LastdbArguments& args,
		                const Alphabet& alph);

		std::istream& readFasta(unsigned numOfIndexes,
		                        const LastdbArguments& args,
		                        const Alphabet& alph,
		                        std::istream& in );

		void accumulateAndFlushPrj();

		void accumulateAndFlushMulti(const LastdbArguments &args);

		void accumulateAndFlushSuffixArrays(int x);

		void createSuffixArrays(int x);

		void replaceVolumeObject();

		void estimateBucketSize();

		static void* threadEntry(void *_args);

//public:
		void startThread();
		void joinThread();

		DatabaseThread(int s, int count);
		~DatabaseThread();
};

void initializeSemaphores();
void destroySemaphores();

#endif //THREADEDLAST_DATABASETHREAD_HH
