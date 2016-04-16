cmake_minimum_required(VERSION 2.6)
project(ThreadedLAST)

set(DIRECTORIES gumbel_params gumbel_params/common
        gumbel_params/corelib)

include_directories(${DIRECTORIES})

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pthread")

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
        prjFiles.hh)

set(LASTAL_SOURCE_FILES
        externalsort.hh
        externalsort.cc
        utilities.hh
        utilities.cc
        heapsort.cc
        heapsort.hh
        linereader.hh
        linereader.cc
        Alphabet.cc
        Alphabet.hh
        MultiSequence.cc
        MultiSequence.hh
        fileMap.hh
        io.hh
        CyclicSubsetSeed.cc
        SubsetSuffixArray.cc
        SubsetSuffixArray.hh
        VectorOrMmap.hh
        LastalArguments.cc
        LastalArguments.hh
        SequenceFormat.hh
        io.cc
        fileMap.cc
        ScoreMatrix.cc
        ScoreMatrix.hh
        DiagonalTable.cc
        DiagonalTable.hh
        SegmentPair.cc
        Alignment.cc
        Alphabet.hh
        GeneralizedAffineGapCosts.hh
        GeneticCode.hh
        GappedXdropAligner.cc
        ScoreMatrixRow.hh
        SegmentPairPot.cc
        SegmentPairPot.hh
        SubsetSuffixArray.cc
        SubsetSuffixArray.hh
        CyclicSubsetSeed.hh
        VectorOrMmap.hh
        Mmap.hh
        fileMap.hh
        AlignmentPot.cc
        AlignmentPot.hh
        Alignment.hh
        GeneralizedAffineGapCosts.cc
        Centroid.cc
        Centroid.hh
        GappedXdropAligner.hh
        GeneralizedAffineGapCosts.hh
        SegmentPair.hh
        GappedXdropAlignerInl.hh
        LambdaCalculator.cc
        LambdaCalculator.hh
        lambda_calculator.hh
        lambda_calculator.cc
        TwoQualityScoreMatrix.cc
        TwoQualityScoreMatrix.hh
        OneQualityScoreMatrix.cc
        OneQualityScoreMatrix.hh
        QualityPssmMaker.cc
        QualityPssmMaker.hh
        qualityScoreUtil.hh
        GeneticCode.cc
        gaplessXdrop.cc
        gaplessXdrop.hh
        io.cc io.hh
        gaplessPssmXdrop.cc
        gaplessPssmXdrop.hh
        gaplessTwoQualityXdrop.cc
        AlignmentWrite.cc
        stringify.hh
        MultiSequenceQual.cc
        MultiSequence.hh
        GappedXdropAlignerPssm.cc
        GappedXdropAligner2qual.cc
        GappedXdropAligner3frame.cc
        lastal.cc
        LastalArguments.hh
        SequenceFormat.hh
        SubsetSuffixArray.hh
        CyclicSubsetSeed.hh
        Mmap.hh
        AlignmentPot.hh
        SegmentPairPot.hh
        ScoreMatrix.hh
        DiagonalTable.hh
        gaplessTwoQualityXdrop.hh
        nrutil.hh
        nrutil.cc
        ludcmp.cc
        lubksb.cc
        LastexArguments.cc
        lastex.cc
        LastexArguments.hh
        version.hh
        utils.hh
        utils.cc
        Alphabet.cc
        Alphabet.hh
        MultiSequence.hh
        MultiSequence.cc
        io.hh
        io.cc
        semaphores.hh
        tempfiles.cc
        tempfiles.hh
        LastEvaluer.hh
        LastEvaluer.cc
        GappedXdropAligner3framePssm.cc
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

add_executable(lastdb ${LASTDB_SOURCE_FILES})
add_executable(lastal ${LASTAL_SOURCE_FILES} ${ALP} ${GUMBEL})