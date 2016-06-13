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

#include "Globals.hh"

using namespace cbrc;

void readPrjFile(
        const std::string &dbname,
        const LastdbArguments &args,
        const Alphabet &alph,
        countT &sequenceCount,
        std::vector<countT> &letterTotals,
        unsigned &volumeNumber);

void generateDifference(const std::string &filename,
                        const std::string &dbname);

void makeAlphabet(Alphabet &alph,
                  const LastdbArguments &args);

bool isDubiousDna(const Alphabet &alph,
                  const MultiSequence &multi);

unsigned makeSubsetSeeds(SubsetSuffixArray indexes[],
                         const LastdbArguments &args,
                         const Alphabet &alph);

void writePrjFile(const std::string &fileName,
                  const LastdbArguments &args,
                  const Alphabet &alph,
                  countT sequenceCount,
                  const std::vector<countT> &letterCounts,
                  unsigned volumes,
                  unsigned numOfIndexes);

void writeOuterPrj(unsigned numOfIndexes);

void renameFiles();

indexT maxLettersPerVolume(const LastdbArguments &args,
                           unsigned numOfIndexes);

void createThreads();

void deleteThreads();

#endif //THREADEDLAST_LASTDB_HH
