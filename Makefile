CXX=g++
CXXFLAGS= -O3 -w -pthread

DBOBJ = src/Alphabet.o \
				src/MultiSequence.o \
				src/CyclicSubsetSeed.o \
				src/SubsetSuffixArray.o \
				src/LastdbArguments.o \
				src/io.o \
				src/fileMap.o  \
				src/SubsetSuffixArraySort.o \
				src/MultiSequenceQual.o \
				src/lastdb.o \
				src/utilities.o \
				src/prjFiles.o \
				src/DatabaseThread.o \
				src/DatabaseWriter.o \
				src/DatabaseVolume.o

ALOBJ = src/lastal.o \
				src/Alphabet.o \
				src/MultiSequence.o \
				src/CyclicSubsetSeed.o \
				src/SubsetSuffixArray.o \
				src/LastalArguments.o \
				src/io.o \
				src/fileMap.o \
				src/ScoreMatrix.o  \
				src/DiagonalTable.o \
				src/SegmentPair.o \
				src/Alignment.o \
				src/GappedXdropAligner.o \
				src/SegmentPairPot.o \
				src/AlignmentPot.o \
				src/GeneralizedAffineGapCosts.o \
				src/Centroid.o \
				src/LambdaCalculator.o \
				src/TwoQualityScoreMatrix.o \
				src/OneQualityScoreMatrix.o \
				src/QualityPssmMaker.o \
				src/GeneticCode.o \
				src/gaplessXdrop.o \
				src/gaplessPssmXdrop.o \
				src/gaplessTwoQualityXdrop.o \
				src/AlignmentWrite.o \
				src/MultiSequenceQual.o \
				src/GappedXdropAlignerPssm.o \
				src/GappedXdropAligner2qual.o \
				src/GappedXdropAligner3frame.o \
				src/lambda_calculator.o \
				src/nrutil.o \
				src/LastexArguments.o \
				src/lastex.o  \
				src/utils.o  \
				src/externalsort.o \
				src/linereader.o \
				src/utilities.o \
				src/heapsort.o \
				src/tempfiles.o \
				src/LastEvaluer.o \
				src/ludcmp.o \
				src/lubksb.o \
				gumbel_params/mcf_local_alignment_evaluer.o \
				gumbel_params/njn_dynprogprob.o \
				gumbel_params/njn_dynprogproblim.o  \
				gumbel_params/njn_dynprogprobproto.o \
				gumbel_params/njn_ioutil.o  \
				gumbel_params/njn_localmaxstat.o \
				gumbel_params/njn_localmaxstatmatrix.o \
				gumbel_params/njn_localmaxstatutil.o \
				gumbel_params/njn_matrix.o \
				gumbel_params/random_gen.o \
				gumbel_params/sls_alp.o \
				gumbel_params/sls_alp_data.o \
				gumbel_params/sls_alp_regression.o \
				gumbel_params/sls_alp_sim.o \
				gumbel_params/sls_pvalues.o \
				alp/sls_alignment_evaluer.o \
				alp/sls_pvalues.o \
				alp/sls_alp_sim.o \
				alp/sls_alp_regression.o \
				alp/sls_alp_data.o \
				alp/sls_alp.o \
				alp/sls_basic.o \
				alp/njn_localmaxstatmatrix.o \
				alp/njn_localmaxstat.o \
				alp/njn_localmaxstatutil.o \
				alp/njn_dynprogprob.o \
				alp/njn_dynprogprobproto.o \
				alp/njn_dynprogproblim.o \
				alp/njn_ioutil.o \
				alp/njn_random.o \
				alp/sls_falp_alignment_evaluer.o \
				alp/sls_fsa1_pvalues.o \
				alp/sls_fsa1_utils.o \
				alp/sls_fsa1.o \
				alp/sls_fsa1_parameters.o

ALL=fastal fastdb

VPATH=src:gumbel_params:alp

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -Igumbel_params -Ialp -Isrc -c -o $@ $<

%.o : %.cc
	$(CXX) $(CXXFLAGS) -Igumbel_params -Ialp -Isrc -c -o $@ $<

all: $(ALL)

fastal: $(ALOBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(ALOBJ)

fastdb: $(DBOBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(DBOBJ) 

clean:
	rm -f */*.o $(ALL)
