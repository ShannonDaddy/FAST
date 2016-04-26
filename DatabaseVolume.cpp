//
// Created by david on 16/04/16.
//

#include "DatabaseVolume.hh"
#include "io.hh"

// Check if we have created a volume that we should proceed to the next volume
/*
bool DatabaseVolume::checkIfReady(){
	return (seqLen > LIMIT);
}
 */

// Take the thread's database structure and add them to the one existing on the disk
/*
void DatabaseVolume::writeToDisk(DatabaseThread *db){

}
 */

void DatabaseVolume::writePooledMultiSequence( const MultiSequence &multi ) const
{
	std::ofstream sspfile ("sspfile", std::ios::binary | std::ios::app );
	std::ofstream tisfile ("tisfile", std::ios::binary | std::ios::app );
	std::ofstream sdsfile ("sdsfile", std::ios::binary | std::ios::app );
	std::ofstream desfile ("desfile", std::ios::binary | std::ios::app );
	std::ofstream quafile ("quafile", std::ios::binary | std::ios::app );

	std::cout << "WRITING THE MULTI FILES TO DISK" << std::endl;
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

	if( multi.qualityScores.begin() != multi.qualityScores.begin() + multi.ends.back() * multi.qualsPerLetter())
	memoryToStream( multi.qualityScores.begin(),
	                multi.qualityScores.begin() + multi.ends.back() * multi.qualsPerLetter(),
	                quafile);
	std::cout << "WRITTEN THE MULTI FILES TO DISK" << std::endl;
}

void DatabaseVolume::writePooledSubsetSuffixArray(const SubsetSuffixArray &sa) const
{
	std::ofstream suffile ("suffile", std::ios::binary | std::ios::app );
	std::ofstream bckfile ("bckfile", std::ios::binary | std::ios::app );

	std::cout << "WRITING THE SUFFIX ARRAY FILES TO DISK" << std::endl;
	if( sa.index.begin() != sa.index.end())
	memoryToStream( sa.index.begin(),
	                sa.index.end(),
	                suffile);
	if( sa.buckets.begin() != sa.buckets.end())
	memoryToStream( sa.buckets.begin(),
	                sa.buckets.end(),
                  bckfile);
	std::cout << "WRITTEN THE SUFFIX ARRAY FILES TO DISK" << std::endl;
}

DatabaseVolume::DatabaseVolume(unsigned _volumes,
                               unsigned _numOfIndexes,
                               unsigned alphSize):
seqLen(0)
{
	prj = new PrjFiles(_volumes, _numOfIndexes, alphSize);
}


DatabaseVolume::~DatabaseVolume(){
	delete prj;
}
