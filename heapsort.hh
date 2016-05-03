#ifndef  __HEAPSORT
#define  __HEAPSORT

#include <vector>
#include <math.h>
#include <utility>
#include "linereader.hh"
using namespace std;

/*
Heap-sort a vector of key-value pairs in descending order of their values.
*/

void heapify(vector< pair<int, Line *> >& A, int i, const std::size_t S);
void build_heap(std::size_t S, vector<pair<int, Line *> >& A);

#endif   // __HEAPSORT
