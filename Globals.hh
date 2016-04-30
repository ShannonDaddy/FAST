//
// Created by david on 25/04/16.
//

#ifndef THREADEDLAST_GLOBALS_HH
#define THREADEDLAST_GLOBALS_HH

#include <string>
#include "Alphabet.hh"
#include "LastdbArguments.hh"
#include "MultiSequence.hh"
#include "DatabaseVolume.hh"
#include <semaphore.h>

using namespace cbrc;

extern Alphabet alph;
extern LastdbArguments args;
extern std::ifstream in;
extern DatabaseVolume *vol;
extern unsigned currentVolumeNumber;

//!! These variables are important
//const std::size_t LOADSIZE = 10000;
const std::size_t LOADSIZE = 100000;
//const unsigned MAXMEMORY = 4294967295; // equivalent to 2^32 - 1 or 4GB.
const unsigned MAXMEMORY = 1048576;
const unsigned maxNumOfIndexes = 16;

typedef MultiSequence::indexT indexT;
typedef unsigned long long countT;

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
#define LOG(x) if( args.verbosity > 0 ) std::cerr << "fastdb: " << x << '\n'

#endif //THREADEDLAST_GLOBALS_HH
