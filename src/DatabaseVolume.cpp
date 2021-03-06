//
// Created by david on 16/04/16.
//

#include "DatabaseVolume.hh"
#include "io.hh"
#include "Globals.hh"

// Check if we have created a volume that we should proceed to the next volume
bool DatabaseVolume::isFinished(unsigned nextBatchEnds,
                                unsigned nextBatchNameEnds) const {

    //std::cout << "endsSize: " << endsSize << std::endl;
    //std::cout << "nameEndsSize: " << nameEndsSize << std::endl;
    //std::cout << (endsSize == nameEndsSize) << std::endl;
    //return endsSize == nameEndsSize;

    return (endsSize + nextBatchEnds < MAXMEMORY && nameEndsSize + nextBatchNameEnds < MAXMEMORY);
}

void DatabaseVolume::writePooledMultiSequence(const MultiSequence &multi,
                                              const LastdbArguments &args,
                                              unsigned endsLastCoordinate,
                                              unsigned nameEndsLastCoordinate) {
    //std::cout << "WRITING THE MULTI FILES TO DISK" << std::endl;
    if (multi.ends.begin() != multi.ends.end())
        memoryToStream(multi.ends.begin() + sizeof(unsigned char),
                       multi.ends.end(),
                       sspfile);

    if (multi.seq.begin() != multi.seq.begin() + endsLastCoordinate)
        memoryToStream(multi.seq.begin() + sizeof(unsigned char),
                       multi.seq.begin() + endsLastCoordinate,
                       tisfile);

    for (int i = 0; i < multi.seq.v.size(); i++) {
        std::cout << multi.seq.v.size() << std::endl;
    }

    if (multi.nameEnds.begin() != multi.nameEnds.begin() + multi.ends.size())
        memoryToStream(multi.nameEnds.begin() + sizeof(unsigned char),
                       multi.nameEnds.begin() + multi.ends.size(),
                       sdsfile);

    if (multi.names.begin() != multi.names.begin() + nameEndsLastCoordinate)
        memoryToStream(multi.names.begin(),
                       multi.names.begin() + nameEndsLastCoordinate,
                       desfile);

    if (args.inputFormat != sequenceFormat::fasta) {
        if (multi.qualityScores.begin() !=
            multi.qualityScores.begin() + multi.ends.back() * multi.qualsPerLetter())
            memoryToStream(multi.qualityScores.begin(),
                           multi.qualityScores.begin() + multi.ends.back() * multi.qualsPerLetter(),
                           quafile);
    }
    //std::cout << "WRITTEN THE MULTI FILES TO DISK" << std::endl;

    endsSize += multi.ends.size();
    nameEndsSize += multi.nameEnds.size();

    endsCoordinate += endsLastCoordinate;
    //!! Is this -1 correct?
    nameEndsCoordinate += nameEndsLastCoordinate - 1;
}

void DatabaseVolume::writePooledSubsetSuffixArray(const SubsetSuffixArray &sa) {
    //std::cout << "WRITING THE SUFFIX ARRAY FILES TO DISK" << std::endl;
    if (sa.index.begin() != sa.index.end())
        memoryToStream(sa.index.begin(),
                       sa.index.end(),
                       suffile);
    //std::cout << "WRITTEN THE SUFFIX ARRAY FILES TO DISK" << std::endl;
}

void DatabaseVolume::writeBucketFile(const SubsetSuffixArray &sa) {
    if (sa.buckets.begin() != sa.buckets.end())
        memoryToStream(sa.buckets.begin(),
                       sa.buckets.end(),
                       bckfile);
}

DatabaseVolume::DatabaseVolume(sequenceFormat::Enum inputFormat,
                               const std::string &dbname,
                               unsigned _volumes,
                               unsigned _numOfIndexes,
                               unsigned alphSize) :
        endsSize(0),
        nameEndsSize(0),
        endsCoordinate(0),
        nameEndsCoordinate(0),
        indexTotal(0) {
    std::stringstream s;
    s << _volumes;
    databaseName = dbname + s.str();
    prj = new PrjFiles(_volumes, _numOfIndexes, alphSize);

    // sequence coordinates (ends)
    // Need to output the first sequence start position (1) to memory in binary format
    sspfile.open((databaseName + ".ssp").c_str(), std::ios::binary);
    int c1 = 1;
    int *a1 = &c1;
    int *b1 = a1 + 1;
    memoryToStream(a1, b1, sspfile);

    // sequences (seq)
    // File starts with an initial pad, assuming padSize = 1
    tisfile.open((databaseName + ".tis").c_str(), std::ios::binary);
    tisfile << ' ';

    // name coordinates (nameEnds)
    // Need to output the first name start position (0) to memory in binary format
    sdsfile.open((databaseName + ".sds").c_str(), std::ios::binary);
    int c2 = 0;
    int *a2 = &c2;
    int *b2 = a2 + 1;
    memoryToStream(a2, b2, sdsfile);

    // names (names)
    desfile.open((databaseName + ".des").c_str(), std::ios::binary);

    if (inputFormat != sequenceFormat::fasta) {
        quafile.open((databaseName + ".qua").c_str(), std::ios::binary);
    }

    suffile.open((databaseName + ".suf").c_str(), std::ios::binary);
    bckfile.open((databaseName + ".bck").c_str(), std::ios::binary);
}

DatabaseVolume::~DatabaseVolume() {
    delete prj;
}

/*
void DatabaseVolume::buildBuckets() {
	indexT bucketDepth = sa.defaultBucketDepth(indexTotal);
	sa.makeBucketSteps( bucketDepth );
	sa.makeBuckets( multi.seqReader(), bucketDepth );
}
 */
