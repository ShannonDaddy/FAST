cmake_minimum_required(VERSION 2.6)
project(FAST)

set(DIRECTORIES
        gumbel_params
        gumbel_params/common
        gumbel_params/corelib
        )

include_directories(${DIRECTORIES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pthread")

set(LASTDB_SOURCE_FILES
        src/Alphabet.cc
        src/Alphabet.hh
        src/MultiSequence.cc
        src/ScoreMatrixRow.hh
        src/CyclicSubsetSeed.cc
        src/SubsetSuffixArray.cc
        src/LastdbArguments.cc
        src/LastdbArguments.hh
        src/SequenceFormat.hh
        src/io.cc
        src/io.hh
        src/fileMap.cc
        src/SubsetSuffixArraySort.cc
        src/SubsetSuffixArraySort.hh
        src/SubsetSuffixArray.hh
        src/CyclicSubsetSeed.hh
        src/Mmap.hh
        src/MultiSequenceQual.cc
        src/MultiSequence.hh
        src/stringify.hh
        src/lastdb.cc
        src/VectorOrMmap.hh
        src/fileMap.hh
        src/qualityScoreUtil.hh
        src/version.hh
        src/lastdb.hh
        src/utilities.cc
        src/utilities.hh
        src/prjFiles.cpp
        src/prjFiles.hh
        src/DatabaseThread.cpp
        src/DatabaseThread.hh
        src/DatabaseWriter.cpp
        src/DatabaseWriter.h
        src/DatabaseVolume.cpp
        src/DatabaseVolume.hh
        src/Globals.hh
        )

set(LASTAL_SOURCE_FILES
        src/externalsort.hh
        src/Alignment.cc
        src/Alignment.hh
        src/AlignmentPot.cc
        src/AlignmentPot.hh
        src/AlignmentWrite.cc
        src/Alphabet.cc
        src/Alphabet.hh
        src/Centroid.cc
        src/Centroid.hh
        src/CyclicSubsetSeed.cc
        src/CyclicSubsetSeed.hh
        src/DiagonalTable.cc
        src/DiagonalTable.hh
        src/externalsort.cc
        src/fileMap.cc
        src/fileMap.hh
        src/gaplessPssmXdrop.cc
        src/gaplessPssmXdrop.hh
        src/gaplessTwoQualityXdrop.cc
        src/gaplessTwoQualityXdrop.hh
        src/gaplessXdrop.cc
        src/gaplessXdrop.hh
        src/GappedXdropAligner.cc
        src/GappedXdropAligner.hh
        src/GappedXdropAligner2qual.cc
        src/GappedXdropAligner3frame.cc
        src/GappedXdropAligner3framePssm.cc
        src/GappedXdropAlignerInl.hh
        src/GappedXdropAlignerPssm.cc
        src/GeneralizedAffineGapCosts.cc
        src/GeneralizedAffineGapCosts.hh
        src/GeneticCode.cc
        src/GeneticCode.hh
        src/heapsort.cc
        src/heapsort.hh
        src/io.cc
        src/io.hh
        src/lambda_calculator.cc
        src/lambda_calculator.hh
        src/LambdaCalculator.cc
        src/LambdaCalculator.hh
        src/lastal.cc
        src/LastalArguments.cc
        src/LastalArguments.hh
        src/LastEvaluer.cc
        src/LastEvaluer.hh
        src/lastex.cc
        src/LastexArguments.cc
        src/LastexArguments.hh
        src/linereader.cc
        src/linereader.hh
        src/lubksb.cc
        src/ludcmp.cc
        src/Mmap.hh
        src/MultiSequence.cc
        src/MultiSequence.hh
        src/MultiSequenceQual.cc
        src/nrutil.cc
        src/nrutil.hh
        src/OneQualityScoreMatrix.cc
        src/OneQualityScoreMatrix.hh
        src/QualityPssmMaker.cc
        src/QualityPssmMaker.hh
        src/qualityScoreUtil.hh
        src/ScoreMatrix.cc
        src/ScoreMatrix.hh
        src/ScoreMatrixRow.hh
        src/SegmentPair.cc
        src/SegmentPair.hh
        src/SegmentPairPot.cc
        src/SegmentPairPot.hh
        src/semaphores.hh
        src/SequenceFormat.hh
        src/stringify.hh
        src/SubsetSuffixArray.cc
        src/SubsetSuffixArray.hh
        src/tempfiles.cc
        src/tempfiles.hh
        src/TwoQualityScoreMatrix.cc
        src/TwoQualityScoreMatrix.hh
        src/utilities.cc
        src/utilities.hh
        src/utils.cc
        src/utils.hh
        src/VectorOrMmap.hh
        src/VectorOrMmap.hh
        src/version.hh
        )

file(GLOB ALP
        "alp/*.h*"
        "alp/*.c*"
        )

file(GLOB GUMBEL
        "gumbel_params/*.h*"
        "gumbel_params/*.c*"
        )

include_directories(${CMAKE_SOURCE_DIR}/alp)
include_directories(${CMAKE_SOURCE_DIR}/gumbel_params)

add_executable(fastdb ${LASTDB_SOURCE_FILES})
add_executable(fastal ${LASTAL_SOURCE_FILES} ${ALP} ${GUMBEL})