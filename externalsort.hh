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

std::vector<std::string> merge_some_files(const std::vector<std::string> &mergelist, 
                                          std::vector<TempFiles> &directories,
                                          const string &tmpdir);

void merge_sorted_files(const vector<string> &filenames,
                        const string &sorted_file_name,
                        const string &tmpdir);

#endif // _EXTERNAL_SORT
