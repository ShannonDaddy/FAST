//
// Created by david on 16/04/16.
//

#ifndef THREADEDLAST_DATABASEVOLUME_HH
#define THREADEDLAST_DATABASEVOLUME_HH

#include <fstream>
#include "prjFiles.hh"
#include "SubsetSuffixArray.hh"

typedef unsigned indexT;
/*
 * File which holds the incomplete database volume until it's ready to be flushed
 */
class DatabaseVolume {
public:
//private:

	std::size_t endsSize;
	std::size_t nameEndsSize;
	indexT endsCoordinate;
	indexT nameEndsCoordinate;

	PrjFiles *prj;

	std::string databaseName;

	std::ofstream sspfile;
	std::ofstream tisfile;
	std::ofstream sdsfile;
	std::ofstream desfile;
	std::ofstream quafile;
	std::ofstream suffile;
	std::ofstream bckfile;

public:
	bool isFinished() const;

	void writePooledMultiSequence( const MultiSequence &multi,
	                               const LastdbArguments &args,
	                               unsigned endsLastCoordinate,
	                               unsigned nameEndsLastCoordinate );

	void writePooledSubsetSuffixArray( const SubsetSuffixArray &sa ) ;

	DatabaseVolume(sequenceFormat::Enum inputFormat,
	               const std::string &dbname,
	               unsigned _volumes,
	               unsigned _numOfIndexes,
	               unsigned alphSize);

	~DatabaseVolume();
};

#endif //THREADEDLAST_DATABASEVOLUME_HH
