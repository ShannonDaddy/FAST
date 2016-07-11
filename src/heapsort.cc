#include "heapsort.hh"

using namespace std;

void heapify(vector<pair<int, Line *> > &A, int i, const std::size_t S) {
    while (true) {
        int l = (2 * i) + 1;
        int r = l + 1;
        int max = i;


        if (l < S && A[l].second->orfid < A[max].second->orfid) {
            max = l;
        } else if (l < S && A[l].second->orfid == A[max].second->orfid) {
            if (l < S && A[l].second->evalue < A[max].second->evalue) {
                max = l;
            } else if (l < S && A[l].second->evalue == A[max].second->evalue) {
                if (l < S && A[l].second->bitscore > A[max].second->bitscore) {
                    max = l;
                }
            }
        }

        if (r < S && A[r].second->orfid < A[max].second->orfid) {
            max = r;
        } else if (r < S && A[r].second->orfid == A[max].second->orfid) {
            if (r < S && A[r].second->evalue < A[max].second->evalue) {
                max = r;
            } else if (r < S && A[r].second->evalue == A[max].second->evalue) {
                if (r < S && A[r].second->bitscore > A[max].second->bitscore) {
                    max = r;
                }
            }
        }

        if (max != i && i < S) {
            pair<int, Line *> temp = A[i];
            A[i] = A[max];
            A[max] = temp;
        } else {
            break;
        }
        i = max;
    }
}

void build_heap(const std::size_t S, vector<pair<int, Line *> > &A) {
    int i = floor(S / 2);
    while (i >= 0) {
        heapify(A, i, S);
        i--;
    }
}