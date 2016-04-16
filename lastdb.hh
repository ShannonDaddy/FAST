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

using namespace cbrc;

#define ERR(x) throw std::runtime_error(x)
#define LOG(x) if( args.verbosity > 0 ) std::cerr << "lastdb: " << x << '\n'

// Defined in lastdb.cc
extern Alphabet alph;
extern LastdbArguments args;
extern std::string currFile;
extern std::ifstream in;

const std::size_t LOADSIZE = 100000;

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



void makeAlphabet( Alphabet& alph, const LastdbArguments& args );
bool isDubiousDna( const Alphabet& alph, const MultiSequence& multi );
unsigned makeSubsetSeeds( SubsetSuffixArray indexes[],
                          const LastdbArguments& args,
                          const Alphabet& alph );

void writePrjFile( const std::string& fileName, const LastdbArguments& args,
                   const Alphabet& alph, countT sequenceCount,
                   const std::vector<countT>& letterCounts,
                   unsigned volumes, unsigned numOfIndexes );

indexT maxLettersPerVolume( const LastdbArguments& args,
                            unsigned numOfIndexes );

#endif //THREADEDLAST_LASTDB_HH
