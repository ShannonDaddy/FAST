//
// Created by david on 16/04/16.
//

#ifndef THREADEDLAST_DATABASEVOLUME_HH
#define THREADEDLAST_DATABASEVOLUME_HH

#include "prjFiles.hh"
#include "SubsetSuffixArray.hh"

/*
 * File which holds the incomplete database volume until it's ready to be flushed
 */
class DatabaseVolume {
public:
//private:
		//MultiSequence multi;
		//SubsetSuffixArray suffixArray;
		PrjFiles *prj;

		unsigned seqLen;

public:
		//bool checkIfReady();
		//void writeToDisk(DatabaseThread *db);
		void writePooledMultiSequence( const MultiSequence &multi ) const;

		void writePooledSubsetSuffixArray(const SubsetSuffixArray &sa) const;

		DatabaseVolume(unsigned _volumes,
		               unsigned _numOfIndexes,
		               unsigned alphSize);
		~DatabaseVolume();
};

#endif //THREADEDLAST_DATABASEVOLUME_HH
