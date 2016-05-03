#ifndef  _EXTERNAL_SORT
#define  _EXTERNAL_SORT

#include <map>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sys/stat.h>

#include "heapsort.hh"
#include "linereader.hh"
#include "utilities.hh"
#include "lastal.hh"
#include "tempfiles.hh"

/*
   All functions required for the initial sequence sorting and blocking
   */

using namespace std;

typedef unsigned long long countT;

void disk_sort_file(const string &outputdir,
                    const string &tobe_sorted_file_name,
                    const string &sorted_file_name,
                    const std::vector<std::string> &mergelist);

std::vector<std::string> mergeFilesInBatches(const std::vector<std::string> &mergelist,
                                             const string &tmpdir);

void merge_sorted_files(const vector<string> &filenames,
                        const string &sorted_file_name,
                        const string &tmpdir);

void openFileHandlers(const vector<string> &filenames,
                      vector<istream_iterator<Line> > &f_its,
                      vector<ifstream *> &ifstream_for_filenames,
                      vector<pair<int, Line *> > &values,
                      Line *curr_lines);

void heapSort(const string &sorted_file_name,
              vector<pair<int, Line *> > &values,
              vector<istream_iterator<Line> > &f_its,
              Line *curr_lines);

void removeFiles(const std::vector<std::string> &batch);

#endif // _EXTERNAL_SORT
