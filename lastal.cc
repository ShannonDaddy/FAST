// Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Martin C. Frith
// BLAST-like pair-wise sequence alignment, using suffix arrays.

#include "lastal.hh"

std::vector<threadData *> *threadDatas;
std::vector<pthread_t> *threads;
std::queue<int> waiting;
std::set<int> working;

SEM_T ioSema;
SEM_T workingSema;
SEM_T waitingQueueSema;
SEM_T workingQueueSema;

unsigned volumes = unsigned(-1);

countT refSequences = -1;
countT refLetters = -1;

void threadData::prepareThreadData(std::string matrixFile, int identifier) {

	output = new outputStruct();

	//unsigned volumes = unsigned(-1);
	indexT minSeedLimit = 0;

	readOuterPrj(args.lastdbName + ".prj", volumes, minSeedLimit, refSequences, refLetters);

	if (minSeedLimit > 1) {
		if (args.outputType == 0)
			ERR("can't use option -j 0: need to re-run lastdb with i <= 1");
		if (minSeedLimit > args.oneHitMultiplicity)
			ERR("can't use option -m < " + stringify(minSeedLimit) +
			    ": need to re-run lastdb with i <= " +
			    stringify(args.oneHitMultiplicity));
	}

	bool isMultiVolume = (volumes + 1 > 0 && volumes > 1);
	args.setDefaultsFromAlphabet(alph.letters == alph.dna, alph.isProtein(),
	                             isCaseSensitiveSeeds, isMultiVolume);
	makeScoreMatrix(matrixFile);

	gapCosts.assign(args.gapExistCost, args.gapExtendCost,
	                args.insExistCost, args.insExtendCost, args.gapPairCost);
	if (args.outputType > 0) {
		calculateScoreStatistics();
	}
	args.setDefaultsFromMatrix(lambdaCalculator.lambda());
	minScoreGapless = args.calcMinScoreGapless(refLetters, numOfIndexes);
	if (!isMultiVolume) {
		args.minScoreGapless = minScoreGapless;
	}
	if (args.outputType > 0) {
		makeQualityScorers();
	}

	if (args.isTranslated()) {
		if (alph.letters == alph.dna) { // allow user-defined alphabet
			ERR("expected protein database, but got DNA");
		}
		queryAlph.fromString(queryAlph.dna);
		if (args.geneticCodeFile.empty()) {
			geneticCode.fromString(geneticCode.standard);
		} else {
			geneticCode.fromFile(args.geneticCodeFile);
		}
		geneticCode.codeTableSet(alph, queryAlph);
		query.initForAppending(3);
	} else {
		queryAlph = alph;
		query.initForAppending(1);
	}
	queryAlph.tr(query.seqWriter(), query.seqWriter() + query.unfinishedSize());

	if (volumes + 1 == 0) {
		readIndex(args.lastdbName, refSequences);
	}

	centroid = new Centroid(gappedXdropAligner);
	this->identifier = identifier;
}

// Set up a scoring matrix, based on the user options
void threadData::makeScoreMatrix(const std::string &matrixFile) {
	if (!matrixFile.empty()) {
		scoreMatrix.fromString(matrixFile);
	}
	else if (args.matchScore < 0 && args.mismatchCost < 0 && alph.isProtein()) {
		const char *b = args.isTranslated() ? "BLOSUM80" : "BLOSUM62";
		scoreMatrix.fromString(ScoreMatrix::stringFromName(b));
	}
	else {
		scoreMatrix.matchMismatch(args.matchScore, args.mismatchCost, alph.letters);
	}

	scoreMatrix.init(alph.encode);

	// If the input is a PSSM, the score matrix is not used, and its
	// maximum score should not be used.  Here, we try to set it to a
	// high enough value that it has no effect.  This is a kludge - it
	// would be nice to use the maximum PSSM score.
	if (args.inputFormat == sequenceFormat::pssm) scoreMatrix.maxScore = 10000;
	// This would work, except the maxDrops aren't finalized yet:
	// maxScore = std::max(args.maxDropGapped, args.maxDropFinal) + 1;
}

void threadData::makeQualityScorers() {
	const ScoreMatrixRow *m = scoreMatrix.caseSensitive;  // case isn't relevant
	double lambda = lambdaCalculator.lambda();
	const std::vector<double> &lp1 = lambdaCalculator.letterProbs1();
	bool isPhred1 = isPhred(referenceFormat);
	int offset1 = qualityOffset(referenceFormat);
	const std::vector<double> &lp2 = lambdaCalculator.letterProbs2();
	bool isPhred2 = isPhred(args.inputFormat);
	int offset2 = qualityOffset(args.inputFormat);

	if (referenceFormat == sequenceFormat::fasta) {
		if (isFastq(args.inputFormat)) {
			if (args.maskLowercase > 0)
				oneQualityScoreMatrixMasked.init(m, alph.size, lambda,
				                                 &lp2[0], isPhred2, offset2,
				                                 alph.canonical, true);
			if (args.maskLowercase < 3)
				oneQualityScoreMatrix.init(m, alph.size, lambda,
				                           &lp2[0], isPhred2, offset2,
				                           alph.canonical, false);
			if (args.outputType > 3) {
				const OneQualityScoreMatrix &m = (args.maskLowercase < 3) ?
				                                 oneQualityScoreMatrix : oneQualityScoreMatrixMasked;
				oneQualityExpMatrix.init(m, args.temperature);
			}
		} else if (args.inputFormat == sequenceFormat::prb) {
			bool isMatchMismatch = (args.matrixFile.empty() && args.matchScore > 0);
			qualityPssmMaker.init(m, alph.size, lambda, isMatchMismatch,
			                      args.matchScore, -args.mismatchCost,
			                      offset2, alph.canonical);
		}
	} else {
		if (isFastq(args.inputFormat)) {
			if (args.maskLowercase > 0)
				twoQualityScoreMatrixMasked.init(m, lambda, &lp1[0], &lp2[0],
				                                 isPhred1, offset1, isPhred2, offset2,
				                                 alph.canonical, true);
			if (args.maskLowercase < 3)
				twoQualityScoreMatrix.init(m, lambda, &lp1[0], &lp2[0],
				                           isPhred1, offset1, isPhred2, offset2,
				                           alph.canonical, false);
			if (args.outputType > 3) {
				ERR("fastq-versus-fastq column probabilities not implemented");
			}
		} else {
			ERR("when the reference is fastq, the query must also be fastq");
		}
	}
}

// Calculate statistical parameters for the alignment scoring scheme
// Meaningless for PSSMs, unless they have the same scale as the score matrix
void threadData::calculateScoreStatistics() {
	LOG("calculating matrix probabilities...");
	// the case-sensitivity of the matrix makes no difference here
	lambdaCalculator.calculate(scoreMatrix.caseSensitive, alph.size);
	if (lambdaCalculator.isBad()) {
		if (isQuality(args.inputFormat) ||
		    (args.temperature < 0 && args.outputType > 3))
			ERR("can't get probabilities for this score matrix");
		else LOG("can't get probabilities for this score matrix");
	} else {
		LOG("lambda=" << lambdaCalculator.lambda());
	}
}

// Read the .prj file for the whole database
void threadData::readOuterPrj(const std::string &fileName, unsigned &volumes, indexT &minSeedLimit,
                              countT &refSequences, countT &refLetters) {

	std::ifstream f(fileName.c_str());
	if (!f) ERR("can't open file: " + fileName);
	unsigned version = 0;

	std::string line, word;
	while (getline(f, line)) {
		std::istringstream iss(line);
		getline(iss, word, '=');
		if (word == "version") iss >> version;
		if (word == "alphabet") iss >> alph;
		if (word == "numofsequences") iss >> refSequences;
		if (word == "numofletters") iss >> refLetters;
		if (word == "maxunsortedinterval") iss >> minSeedLimit;
		if (word == "masklowercase") iss >> isCaseSensitiveSeeds;
		if (word == "sequenceformat") iss >> referenceFormat;
		if (word == "volumes") iss >> volumes;
		if (word == "numofindexes") iss >> numOfIndexes;
	}

	if (f.eof() && !f.bad()) f.clear();
	if (alph.letters.empty() || refSequences + 1 == 0 || refLetters + 1 == 0 ||
	    isCaseSensitiveSeeds < 0 || referenceFormat >= sequenceFormat::prb ||
	    numOfIndexes > maxNumOfIndexes) {
		f.setstate(std::ios::failbit);
	}
	if (!f) ERR("can't read file: " + fileName);
	if (version < 294 && version > 0)
		ERR("the lastdb files are old: please re-run lastdb");
}

// Read a per-volume .prj file, with info about a database volume
void threadData::readInnerPrj(const std::string &fileName,
                              indexT &seqCount, indexT &seqLen) {
	std::ifstream f(fileName.c_str());
	if (!f) ERR("can't open file: " + fileName);

	std::string line, word;
	while (getline(f, line)) {
		std::istringstream iss(line);
		getline(iss, word, '=');
		if (word == "numofsequences") iss >> seqCount;
		if (word == "numofletters") iss >> seqLen;
		if (word == "numofindexes") iss >> numOfIndexes;
	}

	if (f.eof() && !f.bad()) f.clear();
	if (seqCount + 1 == 0 || seqLen + 1 == 0 || numOfIndexes > maxNumOfIndexes) {
		f.setstate(std::ios::failbit);
	}
	if (!f) ERR("can't read file: " + fileName);
}

// Write match counts for each query sequence
void threadData::writeCounts() {

	LOG("writing...");
	std::stringstream outstream;
	std::string output;

	for (indexT i = 0; i < matchCounts.size(); ++i) {
		outstream << query.seqName(i) << "\n";

		for (indexT j = args.minHitDepth; j < matchCounts[i].size(); ++j) {
			outstream << j << "\t" << matchCounts[i][j] << "\n";
		}

		outstream << "\n";  // blank line afterwards
		output = outstream.str();
	}
}

// Count all matches, of all sizes, of a query batch against a suffix array
void threadData::countMatches(char strand) {
	LOG("counting...");
	indexT seqNum = strand == '+' ? 0 : query.finishedSequences() - 1;

	for (indexT i = 0; i < query.finishedSize(); i += args.queryStep) {
		if (strand == '+') {
			for (; ;) {
				if (seqNum == query.finishedSequences()) return;
				if (query.seqEnd(seqNum) > i) break;
				++seqNum;
			}
			// speed-up:
			if (args.minHitDepth > query.seqEnd(seqNum) - i) continue;
		}
		else {
			indexT j = query.finishedSize() - i;
			for (; ;) {
				if (seqNum + 1 == 0) return;
				if (query.seqBeg(seqNum) < j) break;
				--seqNum;
			}
			// speed-up:
			if (args.minHitDepth > j - query.seqBeg(seqNum)) continue;
		}

		for (unsigned x = 0; x < numOfIndexes; ++x)
			subsetUser.countMatches(matchCounts[seqNum], query.seqReader() + i, text.seqReader(),
			                        suffixArrays[x]);
	}
}

// Find query matches to the suffix array, and do gapless extensions
void threadData::alignGapless(SegmentPairPot &gaplessAlns, char strand) {

	Dispatcher dis(Phase::gapless, text, query, scoreMatrix, twoQualityScoreMatrix,
	               twoQualityScoreMatrixMasked,
	               referenceFormat, alph);
	DiagonalTable dt;  // record already-covered positions on each diagonal
	countT matchCount = 0, gaplessExtensionCount = 0, gaplessAlignmentCount = 0;

	for (indexT i = 0; i < query.finishedSize(); i += args.queryStep) {

		for (unsigned x = 0; x < numOfIndexes; ++x) {

			const indexT *beg;
			const indexT *end;

/*
      suffixArrays[x].match( beg, end, dis.b + i, dis.a, args.oneHitMultiplicity, args.minHitDepth );
      matchCount += end - beg;
*/
			subsetUser.match(beg, end, dis.b + i, dis.a, args.oneHitMultiplicity, args.minHitDepth,
			                 suffixArrays[x]);
			matchCount += end - beg;

			// Tried: if we hit a delimiter when using contiguous seeds, then
			// increase "i" to the delimiter position.  This gave a speed-up
			// of only 3%, with 34-nt tags.

			indexT gaplessAlignmentsPerQueryPosition = 0;

			for ( /* no-op*/; beg < end; ++beg) { // loop over suffix-array matches

				if (gaplessAlignmentsPerQueryPosition == args.maxGaplessAlignmentsPerQueryPosition) break;

				indexT j = *beg;  // coordinate in the reference sequence

				if (dt.isCovered(i, j)) continue;

				int fs = dis.forwardGaplessScore(j, i);
				int rs = dis.reverseGaplessScore(j, i);
				int score = fs + rs;
				++gaplessExtensionCount;

				// Tried checking the score after isOptimal & addEndpoint, but
				// the number of extensions decreased by < 10%, and it was
				// slower overall.
				if (score < minScoreGapless) continue;

				indexT tEnd = dis.forwardGaplessEnd(j, i, fs);
				indexT tBeg = dis.reverseGaplessEnd(j, i, rs);
				indexT qBeg = i - (j - tBeg);
				if (!dis.isOptimalGapless(tBeg, tEnd, qBeg)) continue;
				SegmentPair sp(tBeg, qBeg, tEnd - tBeg, score);

				if (args.outputType == 1) {  // we just want gapless alignments
					Alignment aln(identifier);
					aln.fromSegmentPair(sp);
					aln.write(text, query, strand, args.isTranslated(), alph, args.outputFormat, args, output);
				}
				else {
					gaplessAlns.add(sp);  // add the gapless alignment to the pot
				}

				++gaplessAlignmentsPerQueryPosition;
				++gaplessAlignmentCount;
				dt.addEndpoint(sp.end2(), sp.end1());
			}
		}
	}
	LOG("initial matches=" << matchCount);
	LOG("gapless extensions=" << gaplessExtensionCount);
	LOG("gapless alignments=" << gaplessAlignmentCount);
}

// Shrink the SegmentPair to its longest run of identical matches.
// This trims off possibly unreliable parts of the gapless alignment.
// It may not be the best strategy for protein alignment with subset
// seeds: there could be few or no identical matches...
void Dispatcher::shrinkToLongestIdenticalRun(SegmentPair &sp) {
	sp.maxIdenticalRun(a, b, aa->canonical);
	sp.score = gaplessScore(sp.beg1(), sp.end1(), sp.beg2());
}

// Do gapped extensions of the gapless alignments
void threadData::alignGapped(AlignmentPot &gappedAlns, SegmentPairPot &gaplessAlns, Phase::Enum phase) {

	//Dispatcher dis(phase);
	Dispatcher dis(phase, text, query, scoreMatrix, twoQualityScoreMatrix, twoQualityScoreMatrixMasked,
	               referenceFormat, alph);
	indexT frameSize = args.isTranslated() ? (query.finishedSize() / 3) : 0;
	countT gappedExtensionCount = 0, gappedAlignmentCount = 0;

	// Redo the gapless extensions, using gapped score parameters.
	// Without this, if we self-compare a huge sequence, we risk getting
	// huge gapped extensions.
	for (size_t i = 0; i < gaplessAlns.size(); ++i) {
		SegmentPair &sp = gaplessAlns.items[i];

		int fs = dis.forwardGaplessScore(sp.beg1(), sp.beg2());
		int rs = dis.reverseGaplessScore(sp.beg1(), sp.beg2());
		indexT tEnd = dis.forwardGaplessEnd(sp.beg1(), sp.beg2(), fs);
		indexT tBeg = dis.reverseGaplessEnd(sp.beg1(), sp.beg2(), rs);
		indexT qBeg = sp.beg2() - (sp.beg1() - tBeg);
		sp = SegmentPair(tBeg, qBeg, tEnd - tBeg, fs + rs);

		if (!dis.isOptimalGapless(tBeg, tEnd, qBeg)) {
			SegmentPairPot::mark(sp);
		}
	}

	erase_if(gaplessAlns.items, SegmentPairPot::isMarked);

	gaplessAlns.cull(args.cullingLimitForGaplessAlignments);
	gaplessAlns.sort();  // sort by score descending, and remove duplicates

	LOG("redone gapless alignments=" << gaplessAlns.size());

	for (size_t i = 0; i < gaplessAlns.size(); ++i) {
		SegmentPair &sp = gaplessAlns.get(i);

		if (SegmentPairPot::isMarked(sp)) continue;

		Alignment aln(identifier);
		AlignmentExtras extras;  // not used
		aln.seed = sp;

		dis.shrinkToLongestIdenticalRun(aln.seed);

		// do gapped extension from each end of the seed:
		// third...
		aln.makeXdrop(gappedXdropAligner, *centroid, dis.a, dis.b, args.globality,
		              dis.m, scoreMatrix.maxScore, gapCosts, dis.d,
		              args.frameshiftCost, frameSize, dis.p,
		              dis.t, dis.i, dis.j, alph, extras);
		++gappedExtensionCount;

		if (aln.score < args.minScoreGapped) continue;

		if (!aln.isOptimal(dis.a, dis.b, args.globality, dis.m, dis.d, gapCosts,
		                   args.frameshiftCost, frameSize, dis.p,
		                   dis.t, dis.i, dis.j)) {
			// If retained, non-"optimal" alignments can hide "optimal"
			// alignments, e.g. during non-redundantization.
			continue;
		}

		gaplessAlns.markAllOverlaps(aln.blocks);
		gaplessAlns.markTandemRepeats(aln.seed, args.maxRepeatDistance);

		if (phase == Phase::final) gappedAlns.add(aln);
		else SegmentPairPot::markAsGood(sp);

		++gappedAlignmentCount;
	}

	LOG("gapped extensions=" << gappedExtensionCount);
	LOG("gapped alignments=" << gappedAlignmentCount);
}

// Print the gapped alignments, after optionally calculating match
// probabilities and re-aligning using the gamma-centroid algorithm
void threadData::alignFinish(const AlignmentPot &gappedAlns, char strand) {

	Dispatcher dis(Phase::final, text, query, scoreMatrix, twoQualityScoreMatrix,
	               twoQualityScoreMatrixMasked,
	               referenceFormat, alph);
	indexT frameSize = args.isTranslated() ? (query.finishedSize() / 3) : 0;

	if (args.outputType > 3) {
		if (dis.p) {
			LOG("exponentiating PSSM...");
			centroid->setPssm(dis.p, query.finishedSize(), args.temperature,
			                  oneQualityExpMatrix, dis.b, dis.j);
		}
		else {
			centroid->setScoreMatrix(dis.m, args.temperature);
		}
		centroid->setOutputType(args.outputType);
	}

	LOG("finishing...");

	for (size_t i = 0; i < gappedAlns.size(); ++i) {
		const Alignment &aln = gappedAlns.items[i];
		if (args.outputType < 4) {
			aln.write(text, query, strand, args.isTranslated(), alph, args.outputFormat, args, output);
		}
		else {  // calculate match probabilities:
			Alignment probAln(identifier);
			AlignmentExtras extras;
			probAln.seed = aln.seed;
			probAln.makeXdrop(gappedXdropAligner, *centroid,
			                  dis.a, dis.b, args.globality,
			                  dis.m, scoreMatrix.maxScore, gapCosts, dis.d,
			                  args.frameshiftCost, frameSize, dis.p, dis.t,
			                  dis.i, dis.j, alph, extras,
			                  args.gamma, args.outputType);
			probAln.write(text, query, strand, args.isTranslated(),
			              alph, args.outputFormat, args, output, extras);
		}
	}
}

void threadData::makeQualityPssm(bool isApplyMasking) {
	if (!isQuality(args.inputFormat) || isQuality(referenceFormat)) return;

	LOG("making PSSM...");
	query.resizePssm();

	const uchar *seqBeg = query.seqReader();
	const uchar *seqEnd = seqBeg + query.finishedSize();
	const uchar *q = query.qualityReader();
	int *pssm = *query.pssmWriter();

	if (args.inputFormat == sequenceFormat::prb) {
		qualityPssmMaker.make(seqBeg, seqEnd, q, pssm, isApplyMasking);
	}
	else {
		const OneQualityScoreMatrix &m =
				isApplyMasking ? oneQualityScoreMatrixMasked : oneQualityScoreMatrix;
		makePositionSpecificScoreMatrix(m, seqBeg, seqEnd, q, pssm);
	}
}

// Scan one batch of query sequences against one database volume
void threadData::scan(char strand) {

	if (args.outputType == 0) {  // we just want match counts
		countMatches(strand);
		return;
	}

	bool isApplyMasking = (args.maskLowercase > 0);
	makeQualityPssm(isApplyMasking);

	LOG("scanning...");

	SegmentPairPot gaplessAlns;
	alignGapless(gaplessAlns, strand);
	if (args.outputType == 1) return;  // we just want gapless alignments

	if (args.maskLowercase == 1) makeQualityPssm(false);

	AlignmentPot gappedAlns;

	if (args.maskLowercase == 2 || args.maxDropFinal != args.maxDropGapped) {
		alignGapped(gappedAlns, gaplessAlns, Phase::gapped);
		erase_if(gaplessAlns.items, SegmentPairPot::isNotMarkedAsGood);
	}

	if (args.maskLowercase == 2) makeQualityPssm(false);

	alignGapped(gappedAlns, gaplessAlns, Phase::final);

	if (args.outputType > 2) {  // we want non-redundant alignments
		gappedAlns.eraseSuboptimal();
		LOG("nonredundant gapped alignments=" << gappedAlns.size());
	}

	gappedAlns.sort();  // sort by score
	alignFinish(gappedAlns, strand);
}

// Scan one batch of query sequences against one database volume,
// after optionally translating the query
void threadData::translateAndScan(char strand) {

	if (args.isTranslated()) {
		LOG("translating...");
		std::vector<uchar> translation(query.finishedSize());
		geneticCode.translate(query.seqReader(),
		                      query.seqReader() + query.finishedSize(), &translation[0]);

		query.swapSeq(translation);

		scan(strand);
		query.swapSeq(translation);
	}
	else scan(strand);
}

void threadData::readIndex(const std::string &baseName, indexT seqCount) {

	SEM_WAIT(ioSema);
	LOG("reading " << baseName << "...");
	text.fromFiles(baseName, seqCount, isFastq(referenceFormat));
	for (unsigned x = 0; x < numOfIndexes; ++x) {
		if (numOfIndexes > 1) {
			suffixArrays[x].fromFiles(baseName + char('a' + x), isCaseSensitiveSeeds, alph.encode);
		} else {
			suffixArrays[x].fromFiles(baseName, isCaseSensitiveSeeds, alph.encode);
		}
	}
	SEM_POST(ioSema);
}

// Read one database volume
void threadData::readVolume(unsigned volumeNumber) {

	std::string baseName = args.lastdbName + stringify(volumeNumber);
	indexT seqCount = indexT(-1);
	indexT seqLen = indexT(-1);
	readInnerPrj(baseName + ".prj", seqCount, seqLen);
	minScoreGapless = args.calcMinScoreGapless(seqLen, numOfIndexes);
	readIndex(baseName, seqCount);
}

void threadData::reverseComplementPssm() {

	ScoreMatrixRow *beg = query.pssmWriter();
	ScoreMatrixRow *end = beg + query.finishedSize();

	while (beg < end) {
		--end;
		for (unsigned i = 0; i < scoreMatrixRowSize; ++i) {
			unsigned j = queryAlph.complement[i];
			if (beg < end || i < j) std::swap((*beg)[i], (*end)[j]);
		}
		++beg;
	}
}

void threadData::reverseComplementQuery() {
	LOG("reverse complementing...");
	queryAlph.rc(query.seqWriter(), query.seqWriter() + query.finishedSize());
	if (isQuality(args.inputFormat)) {
		std::reverse(query.qualityWriter(),
		             query.qualityWriter() + query.finishedSize() * query.qualsPerLetter());
	} else if (args.inputFormat == sequenceFormat::pssm) {
		reverseComplementPssm();
	}
}

// Scan one batch of query sequences against all database volumes
void threadData::scanAllVolumes(unsigned volumes) {
// Extract this section into the reader function. The threads should just do work as "while buffer non
// empty"
	if (args.outputType == 0) {
		matchCounts.clear();
		matchCounts.resize(query.finishedSequences());
	}

	if (volumes + 1 == 0) volumes = 1;

	for (unsigned i = 0; i < volumes; ++i) {
		if (text.unfinishedSize() == 0 || volumes > 1) readVolume(i);
//
		if (args.strand == 2 && i > 0) reverseComplementQuery();

		if (args.strand != 0) translateAndScan('+');

		if (args.strand == 2 || (args.strand == 0 && i == 0))
			reverseComplementQuery();

		if (args.strand != 1) translateAndScan('-');
	}

	//if( args.outputType == 0 ) writeCounts( out );

	LOG("query batch done!");
}

void writeHeader(countT refSequences, std::ostream &out) {

	out << "# LAST version " <<
	#include "version.hh"
	<< "\n" << "#\n";
	args.writeCommented(out);
	out << "# Reference sequences=" << refSequences << " normal letters=" << refLetters << "\n" << "#\n";

	if (args.outputType == 0) {
		out << "# length\tcount\n";
	} else {
		if (args.outputFormat != 2) {
			if (args.inputFormat != sequenceFormat::pssm || !args.matrixFile.empty()) {
				threadDatas->at(0)->scoreMatrix.writeCommented(out);
				out << "#\n";
			}
			out << "# Coordinates are 0-based.  For - strand matches, coordinates\n";
			out << "# in the reverse complement of the 2nd sequence are used.\n";
			out << "#\n";
		}
		if (args.outputFormat == 0) {  // tabular format
			out << "# score\tname1\tstart1\talnSize1\tstrand1\tseqSize1\t"
			<< "name2\tstart2\talnSize2\tstrand2\tseqSize2\tblocks\n";
		} else if (args.outputFormat == 2) { //blast-like format
			out <<
			"# Q_ID\tS_ID\tIDENT\tALIGN_LEN\tMISMATCHES\tGAPS\tQ_BEG\tQ_END\tS_BEG\tS_END\tE_VAL\tBIT_SCORE\n";
		} else {  // MAF format
			out << "# name start alnSize strand seqSize alignment\n";
		}
	}
	out << "#\n";
}

// Read the next sequence, adding it to the MultiSequence
std::istream &threadData::appendFromFasta(std::istream &in) {

	indexT maxSeqLen = args.batchSize;
	if (maxSeqLen < args.batchSize) maxSeqLen = indexT(-1);
	if (query.finishedSequences() == 0) maxSeqLen = indexT(-1);

	size_t oldUnfinishedSize = query.unfinishedSize();

	/**/ if (args.inputFormat == sequenceFormat::fasta) {
		query.appendFromFasta(in, maxSeqLen);
	} else if (args.inputFormat == sequenceFormat::prb) {
		query.appendFromPrb(in, maxSeqLen, queryAlph.size, queryAlph.decode);
	} else if (args.inputFormat == sequenceFormat::pssm) {
		query.appendFromPssm(in, maxSeqLen, queryAlph.encode,
		                     args.maskLowercase > 1);
	} else {
		query.appendFromFastq(in, maxSeqLen);
	}
	if (!query.isFinished() && query.finishedSequences() == 0) {
		ERR("encountered a sequence that's too long");
	}
	// encode the newly-read sequence
	queryAlph.tr(query.seqWriter() + oldUnfinishedSize,
	             query.seqWriter() + query.unfinishedSize());

	if (isPhred(args.inputFormat))  // assumes one quality code per letter:
		checkQualityCodes(query.qualityReader() + oldUnfinishedSize,
		                  query.qualityReader() + query.unfinishedSize(),
		                  qualityOffset(args.inputFormat));

	return in;
}

void initializeEvalueCalulator(const std::string dbPrjFile, ScoreMatrix &scoreMatrix,
                               std::string dbfilePrj) {
	SequenceStatistics _stats1, _stats2;

	fastaFileSequenceStats(dbfilePrj, &_stats1);
	_stats2 = readStats(dbPrjFile);

	cbrc::LastexArguments args1(args.gapExistCost, args.gapExtendCost);
	initializeEvalueCalculator(args1, scoreMatrix, _stats1, _stats2);

	makeEvaluer();
}

void threadData::callReinit() {

	query.reinitForAppending();
}

void *threadFunction(void *args) {

	struct threadData *data = (struct threadData *) args;
	SEM_WAIT(workingSema);

	SEM_WAIT(workingQueueSema);
	working.insert(data->identifier);
	SEM_POST(workingQueueSema);

	data->scanAllVolumes(volumes);
	data->callReinit();
	/*

					 SEM_WAIT( ioSema );
					 for(int j=0; j < data->output->outputVector->size(); j++){
					 out << data->output->outputVector->at( j );
					 }
					 SEM_POST( ioSema );
				data->output->outputVector->clear();
				*/
	SEM_WAIT(workingQueueSema);
	working.erase(data->identifier);
	SEM_POST(workingQueueSema);

	SEM_WAIT(waitingQueueSema);
	waiting.push(data->identifier);
	SEM_POST(waitingQueueSema);

	SEM_POST(workingSema);
}

void *threadFunctionFinish(void *args) {

	struct threadData *data = (struct threadData *) args;
	SEM_WAIT(workingSema);

	SEM_WAIT(workingQueueSema);
	working.insert(data->identifier);
	SEM_POST(workingQueueSema);

	data->scanAllVolumes(volumes);

	SEM_WAIT(workingQueueSema);
	working.erase(data->identifier);
	SEM_POST(workingQueueSema);

	SEM_WAIT(waitingQueueSema);
	waiting.push(data->identifier);
	SEM_POST(waitingQueueSema);

	SEM_POST(workingSema);
}

void initializeThreads() {

	threadDatas = new std::vector<threadData *>();
	threadDatas->reserve(args.threadNum);

	for (int i = 0; i < args.threadNum; i++) {
		threadData *thread_ptr = new threadData();
		threadDatas->push_back(thread_ptr);
	}

	pthread_t thread;
	threads = new std::vector<pthread_t>(args.threadNum, thread);
}

void initializeSemaphores() {

	//!! Initialize the semaphores
#ifdef __APPLE__
  sem_unlink("/ioSema");
  if ( ( ioSema = sem_open("/ioSema", O_CREAT, 0644, 1)) == SEM_FAILED ) {
    perror("sem_open");
    exit(EXIT_FAILURE);
  }
  sem_unlink("/workingSema");
  if ( ( workingSema = sem_open("/workingSema", O_CREAT, 0644, 1)) == SEM_FAILED ) {
    perror("sem_open");
    exit(EXIT_FAILURE);
  }
  sem_unlink("/waitingQueueSema");
  if ( ( waitingQueueSema = sem_open("/waitingQueueSema", O_CREAT, 0644, 1)) == SEM_FAILED ) {
    perror("sem_open");
    exit(EXIT_FAILURE);
  }
  sem_unlink("/workingQueueSema");
  if ( ( workingQueueSema = sem_open("/workingQueueSema", O_CREAT, 0644, 1)) == SEM_FAILED ) {
    perror("sem_open");
    exit(EXIT_FAILURE);
  }
#elif __linux
	sem_init(&ioSema, 0, 1);
	sem_init(&workingSema, 0, args.threadNum);;
	sem_init(&waitingQueueSema, 0, 1);
	sem_init(&workingQueueSema, 0, 1);
#endif
}

void lastal(int argc, char **argv) {

	args.fromArgs(argc, argv);
	std::string matrixFile;

	if (!args.matrixFile.empty()) {
		matrixFile = ScoreMatrix::stringFromName(args.matrixFile);
		args.fromString(matrixFile);  // read options from the matrix file
		args.fromArgs(argc, argv);  // command line overrides matrix file
	}

	initializeThreads();
	initializeSemaphores();

	std::ofstream outFileStream;
	std::ostream &out = openOut(args.outFile, outFileStream);
	out.precision(3);  // print non-integers more compactly
	countT queryBatchCount = 0;

	char defaultInputName[] = "-";
	char *defaultInput[] = {defaultInputName, 0};
	char **inputBegin = argv + args.inputStart;

	for (int i = 0; i < args.threadNum; i++) {
		threadDatas->at(i)->prepareThreadData(matrixFile, i);
	}

	for (char **i = *inputBegin ? inputBegin : defaultInput; *i; ++i) {
		//writeHeader(refSequences, out);
		std::ifstream inFileStream;
		std::istream &in = openIn(*i, inFileStream);
		initializeEvalueCalulator(args.lastdbName + ".prj", threadDatas->at(0)->scoreMatrix, *inputBegin);

		for (int j = 0; j < args.threadNum; j++) {
			waiting.push(j);
		}

		int io1 = waiting.size();
		int io2 = working.size();
		int id;

		do {
			SEM_WAIT(workingSema);

			if (io1 > 0 && in) {
				SEM_WAIT(waitingQueueSema);
				for (int k = 0; k < waiting.size(); k++) {
					id = waiting.front();
					waiting.pop();
					SEM_POST(waitingQueueSema);

					struct threadData *data = threadDatas->at(id);
					SEM_WAIT(ioSema);
					//data->appendFromFasta( in );
					while (data->appendFromFasta(in)) {
						if (!data->query.isFinished()) {
							pthread_create(&threads->at(id), NULL, threadFunction, (void *) data);
							break;
						}
					}
					SEM_POST(ioSema);
					//pthread_create(&threads->at( id ), NULL, threadFunction, (void*) data);
				}
			}

			SEM_WAIT(waitingQueueSema);
			io1 = waiting.size();
			SEM_POST(waitingQueueSema);
			SEM_WAIT(workingQueueSema);
			io2 = working.size();
			SEM_POST(workingQueueSema);

			SEM_POST(workingSema);

			//} while( ( io1 > 0 || io2 > 0 ) && in );
		} while (in);

		for (int j = 0; j < args.threadNum; j++) {
			pthread_join(threads->at(j), NULL);
			struct threadData *data = threadDatas->at(j);
			pthread_create(&threads->at(id), NULL, threadFunctionFinish, (void *) data);
			pthread_join(threads->at(j), NULL);
		}

	}

	if (!flush(out)) {
		ERR("write error");
	}
}

int main(int argc, char **argv) {

	try {
		lastal(argc, argv);
		return EXIT_SUCCESS;
	} catch (const std::bad_alloc &e) {  // bad_alloc::what() may be unfriendly
		std::cerr << "lastal: out of memory\n";
		std::cout << e.what() << std::endl;
		return EXIT_FAILURE;
	} catch (const std::exception &e) {
		std::cerr << "lastal: " << e.what() << '\n';
		std::cout << e.what() << std::endl;
		return EXIT_FAILURE;
	} catch (int i) {
		return i;
	}
}
