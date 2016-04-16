//
// Created by david on 15/04/16.
//

#ifndef THREADEDLAST_PRJFILES_HH
#define THREADEDLAST_PRJFILES_HH

#include "lastal.hh"
#include "LastdbArguments.hh"

/*
 * This structure is because every thread will create a Prj along with it's database partition.
 * Instead of creating a million mini ones we pour them into a structure and produce it when it has
 * reached a mark of size.
 */

class prjFiles {
private:
		countT sequenceCount;
		unsigned volumes;
		unsigned numOfIndexes;
		std::vector<countT> letterTotals;

public:
		void writePrjFile(const LastdbArguments &args);
		void accumulatePrj(const std::vector<countT>& letterCounts,
		                   unsigned _sequenceCount);
		prjFiles(unsigned _volumes,
		         unsigned numOfIndexes,
		         const Alphabet &alph);
};


#endif //THREADEDLAST_PRJFILES_HH
