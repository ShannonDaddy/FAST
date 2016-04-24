//
// Created by david on 15/04/16.
//

#include <numeric>
#include "prjFiles.hh"

void PrjFiles::accumulatePrj(const std::vector<countT>& letterCounts,
                             unsigned _sequenceCount)
{
	sequenceCount += _sequenceCount;
	for(int i=0; i<letterTotals.size(); i++){
		letterTotals[i] += letterCounts[i];
	}
}

void PrjFiles::writePrjFile(const LastdbArguments &args)
{
	countT letterTotal = std::accumulate( letterTotals.begin(),
	                                      letterTotals.end(), countT(0) );
	//!! placeholder
	std::stringstream s;
	s << volumes;
	std::string fileName = args.lastdbName + s.str() + ".prj";

	std::ofstream f( fileName.c_str() );
	f << "version=" <<
	#include "version.hh"
	<< '\n';
	f << "alphabet=" << alph << '\n';
	f << "numofsequences=" << sequenceCount << '\n';
	f << "numofletters=" << letterTotal << '\n';
	f << "letterfreqs=";
	for( unsigned i = 0; i < letterTotals.size(); ++i ){
		if( i > 0 ) f << ' ';
		f << letterTotals[i];
	}
	f << '\n';

	if( !args.isCountsOnly ){
		f << "maxunsortedinterval=" << args.minSeedLimit << '\n';
		f << "masklowercase=" << args.isCaseSensitive << '\n';
		if( args.inputFormat != sequenceFormat::fasta )
			f << "sequenceformat=" << args.inputFormat << '\n';
		if( volumes+1 > 0 ) f << "volumes=" << volumes << '\n';
		else f << "numofindexes=" << numOfIndexes << '\n';
	}

	if( !f ) ERR( "can't write file: " + fileName );
}

/*
 * prjFiles(-1, numOfIndexes, alph) // this is the master one
 * prjFiles(volumes, numOfIndexes, alph) // per volume
 */
PrjFiles::PrjFiles(unsigned _volumes,
                   unsigned _numOfIndexes,
                   const Alphabet &alph):
volumes(_volumes),
numOfIndexes(_numOfIndexes),
sequenceCount(0)
{
	letterTotals.resize(alph.size);
}
