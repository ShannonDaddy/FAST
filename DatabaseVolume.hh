//
// Created by david on 16/04/16.
//

#ifndef THREADEDLAST_DATABASEVOLUME_HH
#define THREADEDLAST_DATABASEVOLUME_HH

#include <fstream>
#include "prjFiles.hh"
#include "SubsetSuffixArray.hh"

/*
 * File which holds the incomplete database volume until it's ready to be flushed
 */
class DatabaseVolume {
public:
//private:
		//MultiSequence multi;

		std::size_t endsSize;
		std::size_t nameEndsSize;

		//SubsetSuffixArray suffixArray;
		PrjFiles *prj;

		std::string databaseName;

		unsigned seqLen;

		std::ofstream sspfile;
		std::ofstream tisfile;
		std::ofstream sdsfile;
		std::ofstream desfile;
		std::ofstream quafile;
		std::ofstream suffile;
		std::ofstream bckfile;

public:
		bool isFinished() const;
		//void writeToDisk(DatabaseThread *db);
		void writePooledMultiSequence( const MultiSequence &multi,
		                               const LastdbArguments &args) ;

		void writePooledSubsetSuffixArray( const SubsetSuffixArray &sa ) ;

		DatabaseVolume(sequenceFormat::Enum inputFormat,
		               const std::string &dbname,
		               unsigned _volumes,
		               unsigned _numOfIndexes,
		               unsigned alphSize);

		~DatabaseVolume();
};

#endif //THREADEDLAST_DATABASEVOLUME_HH
