//
// Created by david on 13/04/16.
//

#ifndef THREADEDLAST_LASTDB_HH
#define THREADEDLAST_LASTDB_HH

#include "LastdbArguments.hh"
#include "SubsetSuffixArray.hh"
#include "Alphabet.hh"
#include "MultiSequence.hh"
#include "io.hh"
#include "qualityScoreUtil.hh"
#include <numeric>  // accumulate

#include <semaphore.h>

#ifdef __APPLE__
typedef sem_t* SEM_T;
  #define SEM_POST(x) sem_post(x)
  #define SEM_WAIT(x) sem_wait(x)
#elif __linux
typedef sem_t SEM_T;
#define SEM_POST(x) sem_post(&x)
#define SEM_WAIT(x) sem_wait(&x)
#endif

#define ERR(x) throw std::runtime_error(x)
#define LOG(x) if( args.verbosity > 0 ) std::cerr << "lastdb: " << x << '\n'

const std::size_t LOADSIZE = 100000;

using namespace cbrc;

typedef MultiSequence::indexT indexT;
typedef unsigned long long countT;

const unsigned maxNumOfIndexes = 16;

void readPrjFile(
		const std::string &dbname,
		const LastdbArguments &args,
		const Alphabet &alph,
		countT &sequenceCount,
		std::vector<countT> &letterTotals,
		unsigned &volumeNumber);

void generateDifference(const std::string &filename, const std::string &dbname);

void initializeSemaphores();
void destroySemaphores();

class DatabaseThread{
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

struct File{

		//std::ifstream inFileStream;
		//std::istream &in = openIn( inputName, inFileStream );
};

#endif //THREADEDLAST_LASTDB_HH
