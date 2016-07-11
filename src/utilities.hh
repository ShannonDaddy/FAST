#ifndef ___UTILITIES___
#define ___UTILITIES___

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <regex.h>
#include <list>

using namespace std;

char *split_n_pick(const string &strn, char *buf, char d, unsigned n);

std::string random_str(int len);

void topHits(const std::string &filename, int maxHits);

//void topHitsVector(std::vector<std::string> &outputVector, int maxHits);
void topHitsVector(std::list<std::string> &outputVector, int maxHits);

string orf_extractor_from_blast(const string &line);

double evalue_extractor_from_blast(const string &line);

double bit_score_extractor_from_blast(const string &line);

string generate_directory_name(const string &tmpdir);

#endif //_UTILITIES
