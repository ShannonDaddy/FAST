#include <iostream>
#include <stack>
#include <stdlib.h>
#include <dirent.h>
#include <fstream>
#include <sys/stat.h>
#include <vector>
#include <cstring>

using namespace std;

class TempFiles {

private:
	string filename(unsigned int i);
	string dirname(unsigned int i);
	string toString(unsigned int i);
	void remove_dir(char *path);
	void setFanOut(unsigned int i ) {  S = i; }

	vector<string> filenames;
	std::string tempdir;
	std::string basedir ;
	unsigned int count;
	unsigned int S; // fanout

public:
	TempFiles(const std::string &_tempdir,
	          const std::string &_basedir);

	string nextFileName() ;
	void clear();
	std::size_t size();
	vector<string> getFileNames();
};
