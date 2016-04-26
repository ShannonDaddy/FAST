//
// Created by david on 16/04/16.
//

#include "DatabaseVolume.hh"

// Check if we have created a volume that we should proceed to the next volume
bool DatabaseVolume::checkIfReady(){
	return (seqLen > LIMIT);
}

// Take the thread's database structure and add them to the one existing on the disk
void DatabaseVolume::writeToDisk(DatabaseThread *db){

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
