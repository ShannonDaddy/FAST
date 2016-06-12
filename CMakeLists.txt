cmake_minimum_required(VERSION 2.6)
project(ThreadedLAST)

set(DIRECTORIES gumbel_params gumbel_params/common
        gumbel_params/corelib)

include_directories(${DIRECTORIES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pthread")

set(LASTDB_SOURCE_FILES
        Alphabet.cc
        Alphabet.hh
        MultiSequence.cc
        ScoreMatrixRow.hh
        CyclicSubsetSeed.cc
        SubsetSuffixArray.cc
        LastdbArguments.cc
        LastdbArguments.hh
        SequenceFormat.hh
        io.cc
        io.hh
        fileMap.cc
        SubsetSuffixArraySort.cc
        SubsetSuffixArraySort.hh
        SubsetSuffixArray.hh
        CyclicSubsetSeed.hh
        Mmap.hh
        MultiSequenceQual.cc
        MultiSequence.hh
        stringify.hh
        lastdb.cc
        VectorOrMmap.hh
        fileMap.hh
        qualityScoreUtil.hh
        version.hh
        lastdb.hh
        utilities.cc
        utilities.hh
        prjFiles.cpp
        prjFiles.hh
        DatabaseThread.cpp
        DatabaseThread.hh
        DatabaseWriter.cpp
        DatabaseWriter.h
        DatabaseVolume.cpp
        DatabaseVolume.hh
        Globals.hh
        )

set(LASTAL_SOURCE_FILES
        externalsort.hh
        Alignment.cc
        Alignment.hh
        AlignmentPot.cc
        AlignmentPot.hh
        AlignmentWrite.cc
        Alphabet.cc
        Alphabet.hh
        Centroid.cc
        Centroid.hh
        CyclicSubsetSeed.cc
        CyclicSubsetSeed.hh
        DiagonalTable.cc
        DiagonalTable.hh
        externalsort.cc
        fileMap.cc
        fileMap.hh
        gaplessPssmXdrop.cc
        gaplessPssmXdrop.hh
        gaplessTwoQualityXdrop.cc
        gaplessTwoQualityXdrop.hh
        gaplessXdrop.cc
        gaplessXdrop.hh
        GappedXdropAligner.cc
        GappedXdropAligner.hh
        GappedXdropAligner2qual.cc
        GappedXdropAligner3frame.cc
        GappedXdropAligner3framePssm.cc
        GappedXdropAlignerInl.hh
        GappedXdropAlignerPssm.cc
        GeneralizedAffineGapCosts.cc
        GeneralizedAffineGapCosts.hh
        GeneticCode.cc
        GeneticCode.hh
        heapsort.cc
        heapsort.hh
        io.cc
        io.hh
        lambda_calculator.cc
        lambda_calculator.hh
        LambdaCalculator.cc
        LambdaCalculator.hh
        lastal.cc
        LastalArguments.cc
        LastalArguments.hh
        LastEvaluer.cc
        LastEvaluer.hh
        lastex.cc
        LastexArguments.cc
        LastexArguments.hh
        linereader.cc
        linereader.hh
        lubksb.cc
        ludcmp.cc
        Mmap.hh
        MultiSequence.cc
        MultiSequence.hh
        MultiSequenceQual.cc
        nrutil.cc
        nrutil.hh
        OneQualityScoreMatrix.cc
        OneQualityScoreMatrix.hh
        QualityPssmMaker.cc
        QualityPssmMaker.hh
        qualityScoreUtil.hh
        ScoreMatrix.cc
        ScoreMatrix.hh
        ScoreMatrixRow.hh
        SegmentPair.cc
        SegmentPair.hh
        SegmentPairPot.cc
        SegmentPairPot.hh
        semaphores.hh
        SequenceFormat.hh
        stringify.hh
        SubsetSuffixArray.cc
        SubsetSuffixArray.hh
        tempfiles.cc
        tempfiles.hh
        TwoQualityScoreMatrix.cc
        TwoQualityScoreMatrix.hh
        utilities.cc
        utilities.hh
        utils.cc
        utils.hh
        VectorOrMmap.hh
        VectorOrMmap.hh
        version.hh
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