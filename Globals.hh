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
extern std::string currFile;
extern std::ifstream in;
extern DatabaseVolume *vol;

const std::size_t LOADSIZE = 100000;

typedef MultiSequence::indexT indexT;
typedef unsigned long long countT;

//!! Absolutely arbitrary number that I need as a placeholder for now
const unsigned LIMIT = 10000000;
const unsigned maxNumOfIndexes = 16;

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
