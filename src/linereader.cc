// FastaLine.cc                                                   -*-C++-*-

#include "linereader.hh"

#include <iomanip>

// CREATORS
Line::Line()
        : orfid(), evalue(), bitscore(), line() {
}

Line::Line(const std::string &name, const std::string &sequence)
        : orfid(name), evalue(evalue_extractor_from_blast(sequence)),
          bitscore(bit_score_extractor_from_blast(sequence)), line(sequence) {
}

Line::Line(const Line &original)
        : orfid(original.orfid), evalue(original.evalue), bitscore(original.bitscore), line(original.line) {
}

// MANIPULATORS
Line &Line::operator=(const Line &rhs) {
    if (this != &rhs) {
        orfid = rhs.orfid;
        evalue = rhs.evalue;
        bitscore = rhs.bitscore;
        line = rhs.line;
    }

    return *this;
}

void Line::setOrfId(const std::string &value) {
    orfid.assign(value.begin(), value.end());
}

void Line::setLine(const std::string &value) {
    line.assign(value.begin(), value.end());
}

void Line::setEvalue(double value) {
    evalue = value;
}

void Line::setBitscore(double value) {
    bitscore = value;
}

// ACCESSORS
const std::string &Line::getOrfId() const {
    return orfid;
}

const std::string &Line::getLine() const {
    return line;
}


// FREE OPERATORS
bool operator==(const Line &lhs, const Line &rhs) {
    return lhs.getOrfId() == rhs.getOrfId() && lhs.getLine() == rhs.getLine();
}

bool operator!=(const Line &lhs, const Line &rhs) {
    return lhs.getOrfId() != rhs.getOrfId() || lhs.getLine() != rhs.getLine();
}

std::ostream &Line::print(std::ostream &stream) const {
    if (stream.good()) {
        stream << this->line << std::endl;
    }

    return stream;
}

// FREE OPERATORS
std::istream &operator>>(std::istream &stream, Line &rhs) {
    char buf[10000];
    std::string _line;
    if (getline(stream, _line)) {
        rhs.line = _line;
        char *field = split_n_pick(_line, buf, '\t', 0);
        rhs.orfid = string(field);
        rhs.evalue = evalue_extractor_from_blast(_line);
        rhs.bitscore = bit_score_extractor_from_blast(_line);
    }
    return stream;
}

std::ostream &operator<<(std::ostream &stream, const Line &rhs) {
    stream << rhs.line;
    return stream;
}

// close package namespace
