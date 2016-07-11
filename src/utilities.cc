#include <sys/stat.h>
#include "utilities.hh"

char *split_n_pick(const string &strn, char *buf, char d, unsigned n) {
    strcpy(buf, strn.c_str());

    char *v = buf;
    char *s1 = buf;
    v = s1;

    unsigned int i = 0;

    while (*s1 != '\0') {
        if (*s1 == d) {
            *s1 = '\0';
            i++;
            if (i > n) return v;
            v = s1 + 1;
        }
        s1++;
    }
    return v;
}

string random_str(int len) {
    static const char alphanum[] =
            "0123456789"
                    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                    "abcdefghijklmnopqrstuvwxyz";

    string str;

    for (int i = 0; i < len; ++i) {
        str += alphanum[rand() % (sizeof(alphanum) - 1)];
    }
    return str;
}

string generate_directory_name(const string &tmpdir) {
    string random_string;
    int err = -1;
    do {
        random_string = random_str(30);
        string potential_directory = tmpdir + random_string + "LASTtemp0";
        struct stat potential_directory_stat;
        err = stat(potential_directory.c_str(), &potential_directory_stat);
    } while (err != -1);
    return random_string;
}

string orf_extractor_from_blast(const string &line) {
    char buf[10000];
    string orfid = split_n_pick(line, buf, '\t', 0);
    return orfid;
}

double evalue_extractor_from_blast(const string &line) {
    char buf[10000];
    string evaluestr = split_n_pick(line, buf, '\t', 10);
    double evalue;

    try {
        evalue = atof(evaluestr.c_str());
    }
    catch (...) {
        return 100;
    }
    return evalue;
}

double bit_score_extractor_from_blast(const string &line) {
    char buf[10000];
    string value = split_n_pick(line, buf, '\t', 11);
    double bitscore;

    try {
        bitscore = atof(value.c_str());
    }
    catch (...) {
        return 0;
    }
    return bitscore;
}


void topHits(const std::string &filename, int maxHits) {
    int count = 0;
    std::ifstream input(filename.c_str());
    std::string current;
    std::string prevorfid = "";
    std::ofstream output((filename + "_tmp").c_str());

    while (getline(input, current)) {
        std::size_t location = current.find_first_of("\t");

        if (!(prevorfid.compare(0, location, current) == 0 || prevorfid.size() == 0))
            count = 0;

        if (count < maxHits) {
            output << current << "\n";
        }

        count++;
        prevorfid.assign(current, 0, location);
    }
    std::rename((filename + "_tmp").c_str(), filename.c_str());
}

void topHitsVector(std::list<std::string> &outputVector, int maxHits) {
    int count = 0;
    std::string prevorfid = "";

    std::list<std::string>::iterator it = outputVector.begin();
    for (; it != outputVector.end(); ++it) {
        std::size_t location = it->find_first_of("\t");

        if (!(it->compare(0, location, *it) == 0 || prevorfid.size() == 0))
            count = 0;

        if (count >= maxHits)
            *it = "";

        count++;
        prevorfid.assign(*it, 0, location);
    }
}
