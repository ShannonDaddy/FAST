//
// Created by david on 16/04/16.
//

#ifndef THREADEDLAST_DATABASEVOLUME_HH
#define THREADEDLAST_DATABASEVOLUME_HH

#include "DatabaseThread.hh"
//!! Recursive include is causing headaches
//#include "prjFiles.hh"

//!! Absolutely arbitrary number that I need as a placeholder for now
const unsigned LIMIT = 10000000;
/*
 * File which holds the incomplete database volume until it's ready to be flushed
 */
class DatabaseVolume {
private:
		MultiSequence multi;
		SubsetSuffixArray suffixArray;
		//PrjFiles prj;

		unsigned seqLen;

public:
		bool checkIfReady();
		void writeToDisk(DatabaseThread *db);

		DatabaseVolume();
};


#endif //THREADEDLAST_DATABASEVOLUME_HH
