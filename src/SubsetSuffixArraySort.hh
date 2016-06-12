//
// Created by david on 14/04/16.
//
#ifndef THREADEDLAST_SUBSETSUFFIXARRAYSORT_HH
#define THREADEDLAST_SUBSETSUFFIXARRAYSORT_HH

#include "SubsetSuffixArray.hh"
#include <algorithm>  // iter_swap, min

using namespace cbrc;

typedef SubsetSuffixArray::indexT indexT;
struct Stack {
		indexT *beg;
		indexT *end;
		indexT depth;
};

namespace cbrc {
	class SuffixArraySorter {

	private:
			Stack stack[1048576];  // big enough???
			Stack *sp;

			friend class SubsetSuffixArray;

			inline void PUSH(indexT *b, indexT *e, indexT d) {
				sp->beg = b;
				sp->end = e;
				(sp++)->depth = d;
			}

			inline void POP(indexT *&b, indexT *&e, indexT &d) {
				b = (--sp)->beg;
				e = sp->end;
				d = sp->depth;
			}

	public:
			void insertionSort(const uchar *text, const CyclicSubsetSeed &seed,
			                   indexT *beg, indexT *end, indexT depth);

			void radixSort1(const uchar *text, const uchar *subsetMap,
			                indexT *beg, indexT *end, indexT depth);

			void radixSort2(const uchar *text, const uchar *subsetMap,
			                indexT *beg, indexT *end, indexT depth);

			void radixSort3(const uchar *text, const uchar *subsetMap,
			                indexT *beg, indexT *end, indexT depth);

			// Specialized sort for 4 symbols + 1 delimiter.  E.g. DNA.
			void radixSort4(const uchar *text, const uchar *subsetMap,
			                indexT *beg, indexT *end, indexT depth);

			void radixSortN(const uchar *text, const uchar *subsetMap,
			                indexT *beg, indexT *end, indexT depth,
			                unsigned subsetCount);

			SuffixArraySorter() { sp = stack; }
	};
}

#endif //THREADEDLAST_SUBSETSUFFIXARRAYSORT_HH
