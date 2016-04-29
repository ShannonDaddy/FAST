//
// Created by david on 16/04/16.
//

#include "DatabaseVolume.hh"
#include "io.hh"

// Check if we have created a volume that we should proceed to the next volume
bool DatabaseVolume::isFinished() const {

	//std::cout << "endsSize: " << endsSize << std::endl;
	//std::cout << "nameEndsSize: " << nameEndsSize << std::endl;
	//std::cout << (endsSize == nameEndsSize) << std::endl;
	return endsSize == nameEndsSize;

	//return true;
}

void DatabaseVolume::writePooledMultiSequence( const MultiSequence &multi,
                                               const LastdbArguments &args)
{

	//std::cout << "WRITING THE MULTI FILES TO DISK" << std::endl;
	if( multi.ends.begin() != multi.ends.end() )
	memoryToStream( multi.ends.begin() + sizeof(unsigned char),
	                multi.ends.end(),
	                sspfile );

	if( multi.seq.begin() != multi.seq.begin() + multi.ends.back())
	memoryToStream( multi.seq.begin() + sizeof(unsigned char),
	                multi.seq.begin() + multi.ends.back(),
	                tisfile );

	if( multi.nameEnds.begin() != multi.nameEnds.begin() + multi.ends.size())
	memoryToStream( multi.nameEnds.begin() + sizeof(unsigned char),
	                multi.nameEnds.begin() + multi.ends.size(),
	                sdsfile);

	if( multi.names.begin() != multi.names.begin() + multi.nameEnds[multi.finishedSequences()] )
	memoryToStream( multi.names.begin(),
	                multi.names.begin() + multi.nameEnds[multi.finishedSequences()],
	                desfile);

	if ( args.inputFormat != sequenceFormat::fasta ) {
		if (multi.qualityScores.begin() !=
		    multi.qualityScores.begin() + multi.ends.back() * multi.qualsPerLetter())
			memoryToStream(multi.qualityScores.begin(),
			               multi.qualityScores.begin() + multi.ends.back() * multi.qualsPerLetter(),
			               quafile);
	}
	//std::cout << "WRITTEN THE MULTI FILES TO DISK" << std::endl;

	endsSize += multi.ends.size();
	nameEndsSize += multi.nameEnds.size();

	endsCoordinate += multi.ends.back();
	nameEndsCoordinate += multi.nameEnds.back() - 1;
}

void DatabaseVolume::writePooledSubsetSuffixArray(const SubsetSuffixArray &sa)
{
	//std::cout << "WRITING THE SUFFIX ARRAY FILES TO DISK" << std::endl;
	if( sa.index.begin() != sa.index.end())
	memoryToStream( sa.index.begin(),
	                sa.index.end(),
	                suffile);
	if( sa.buckets.begin() != sa.buckets.end())
	memoryToStream( sa.buckets.begin(),
	                sa.buckets.end(),
                  bckfile);
	//std::cout << "WRITTEN THE SUFFIX ARRAY FILES TO DISK" << std::endl;
}

DatabaseVolume::DatabaseVolume(sequenceFormat::Enum inputFormat,
                               const std::string &dbname,
                               unsigned _volumes,
                               unsigned _numOfIndexes,
                               unsigned alphSize):
		endsSize(0),
		nameEndsSize(0),
		endsCoordinate(0),
		nameEndsCoordinate(0)
{
	std::stringstream s;
	s << _volumes;
	databaseName = dbname + s.str();
	prj = new PrjFiles(_volumes, _numOfIndexes, alphSize);

	// sequence coordinates (ends)
	// Need to output the first sequence start position (1) to memory in binary format
	sspfile.open( (databaseName+".ssp").c_str(), std::ios::binary );
	int c1 = 1;
	int *a1 = &c1;
	int *b1 = a1+1;
	memoryToStream(a1, b1, sspfile);

	// sequences (seq)
	// File starts with an initial pad, assuming padSize = 1
	tisfile.open( (databaseName+".tis").c_str(), std::ios::binary );
	tisfile << ' ';

	// name coordinates (nameEnds)
	// Need to output the first name start position (0) to memory in binary format
	sdsfile.open( (databaseName+".sds").c_str(), std::ios::binary );
	int c2 = 0;
	int *a2 = &c2;
	int *b2 = a2+1;
	memoryToStream(a2, b2, sdsfile);

	// names (names)
	desfile.open( (databaseName+".des").c_str(), std::ios::binary );

	if ( inputFormat != sequenceFormat::fasta ) {
		quafile.open( (databaseName+".qua").c_str(), std::ios::binary );
	}

	suffile.open( (databaseName+".suf").c_str(), std::ios::binary );
	bckfile.open( (databaseName+".bck").c_str(), std::ios::binary );
}

DatabaseVolume::~DatabaseVolume(){
	delete prj;
}
