//
// Created by david on 16/04/16.
//

#include "DatabaseVolume.hh"
#include "io.hh"

// Check if we have created a volume that we should proceed to the next volume
bool DatabaseVolume::isFinished() const {
	//std::cout << "endsSize: " << endsSize << std::endl;
	//std::cout << "nameEndsSize: " << nameEndsSize << std::endl;
	return endsSize == nameEndsSize;
}

void DatabaseVolume::writePooledMultiSequence( const MultiSequence &multi,
                                               const LastdbArguments &args)
{
	endsSize += multi.ends.size();
	nameEndsSize += multi.nameEnds.size();

	//std::cout << "WRITING THE MULTI FILES TO DISK" << std::endl;
	if( multi.ends.begin() != multi.ends.end() )
	memoryToStream( multi.ends.begin(),
	                multi.ends.end(),
	                sspfile );

	if( multi.seq.begin() != multi.seq.begin() + multi.ends.back())
	memoryToStream( multi.seq.begin(),
	                multi.seq.begin() + multi.ends.back(),
	                tisfile );

	if( multi.nameEnds.begin() != multi.nameEnds.begin() + multi.ends.size())
	memoryToStream( multi.nameEnds.begin(),
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
		seqLen(0),
		endsSize(0),
		nameEndsSize(0)
{
	std::stringstream s;
	s << _volumes;
	databaseName = dbname + s.str();
	prj = new PrjFiles(_volumes, _numOfIndexes, alphSize);

	sspfile.open( (databaseName+".ssp").c_str(), std::ios::binary );
	tisfile.open( (databaseName+".tis").c_str(), std::ios::binary );
	sdsfile.open( (databaseName+".sds").c_str(), std::ios::binary );
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
