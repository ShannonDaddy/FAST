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

		unsigned volumeNumber;
		unsigned numOfIndexes;
		countT sequenceCount;
		std::vector<countT> letterCounts;
		std::vector<countT> letterTotals;

		int rank;
		SuffixArraySorter *sorter;

		void formatdb(const LastdbArguments &args,
		              const Alphabet &alph,
		              const unsigned numOfIndexes,
		              const std::string &inputName);
		void prepareNextVolume();

		void makeVolume( unsigned numOfIndexes,
		                 const LastdbArguments& args,
		                 const Alphabet& alph,
		                 const std::vector<countT>& letterCounts,
		                 const std::string& baseName );

		std::istream& readFasta(unsigned numOfIndexes,
		                        const LastdbArguments& args,
		                        const Alphabet& alph,
		                        std::istream& in );

		static void* threadEntry(void *args);
		void threadFunction();

//public:
		void startThread();
		void joinThread();

		DatabaseThread(int s, int count);
		~DatabaseThread();
};

void initializeSemaphores();
void destroySemaphores();

#endif //THREADEDLAST_DATABASETHREAD_HH
