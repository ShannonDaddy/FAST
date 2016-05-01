//
// Created by david on 15/04/16.
//

#ifndef THREADEDLAST_PRJFILES_HH
#define THREADEDLAST_PRJFILES_HH

#include "LastdbArguments.hh"
#include "MultiSequence.hh"
#include "Alphabet.hh"
#include "SubsetSuffixArray.hh"

/*
 * This structure is because every thread will create a Prj along with it's database partition.
 * Instead of creating a million mini ones we pour them into a structure and produce it when it has
 * reached a mark of size.
 */

using namespace cbrc;

typedef MultiSequence::indexT indexT;
typedef unsigned long long countT;

class PrjFiles {
//private:
public:
	countT sequenceTotal;
	unsigned volumes;
	unsigned numOfIndexes;
	std::vector<countT> letterTotals;

public:
	void writePrjFile(const LastdbArguments &args,
	                  const SubsetSuffixArray &sa,
	                  indexT textLength,
	                  std::size_t indexTotal);

	void accumulatePrj(const std::vector<countT>& letterCounts,
	                   unsigned _sequenceCount);

	PrjFiles(unsigned _volumes,
	         unsigned numOfIndexes,
	         unsigned alphSize);
};


#endif //THREADEDLAST_PRJFILES_HH
