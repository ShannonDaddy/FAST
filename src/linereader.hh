#ifndef INCLUDED_FASTA_Line
#define INCLUDED_FASTA_Line

#include <string>
#include <iterator>
#include "utilities.hh"

using namespace std;

class Line {
    friend std::istream &operator>>(std::istream &, Line &);

public:
    std::string orfid;      // name of the Line
    std::string line;  // sequence of 'A', 'C', 'G', and 'T' characters
    double evalue;
    double bitscore;

    explicit Line();

    Line(const std::string &orfid, const std::string &line);
    // Create a new 'Line' object having the specified 'name' and
    // 'sequence' attribute values.

    Line(const Line &original);
    // Create a new 'Line' object having the same value as the specified
    // 'original' object.

    // MANIPULATORS
    Line &operator=(const Line &rhs);
    // Set the value of this object to the value of the specified 'rhs'
    // object, and return a reference providing modifiable access to this
    // object.

    void setOrfId(const std::string &value);
    // Set the 'name' attribute of this object to the specified 'value'.

    void setLine(const std::string &value);
    // Set the 'sequence' attribute of this object to the specified
    // 'value'.

    void setEvalue(double value);

    void setBitscore(double value);

    // ACCESSORS
    const std::string &getOrfId() const;
    // Return a reference providing non-modifiable access to the 'name'
    // attribute of this object.

    const std::string &getLine() const;
    // Return a reference providing non-modifiable access to the 'sequence'
    // attribute of this object.

    // Aspects

    /*std::ostream& print(std::ostream& stream,
                        int           level = 0,
                        int           spacesPerLevel = 4) const;*/

    std::ostream &print(std::ostream &stream) const;

};

// FREE OPERATORS
bool operator==(const Line &lhs, const Line &rhs);
// Return 'true' if the specified 'lhs' and 'rhs' objects have the same
// value and 'false' otherwise.  Two 'Line' objects have the same value
// if the corresponding values of their 'sequence' and 'name' attributes have
// the same value.

bool operator!=(const Line &lhs, const Line &rhs);
// Return 'true' if the specified 'lhs' and 'rhs' objects do not have the
// same value and 'false' otherwise.  Two 'Line' objects do not have the
// same value if any the corresponding values of their 'sequence' and 'name'
// attributes do not have the same value.

std::istream &operator>>(std::istream &stream, Line &rhs);
// Assign to the specified 'rhs' object the value extracted from the
// specified 'stream', and return a reference providing modifiable access
// to 'stream'.

std::ostream &operator<<(std::ostream &stream, const Line &rhs);
// Output the value of the specified 'rhs' object to the specified
// 'stream', and return a reference providing modifiable access to
// 'stream'.

// close package namespace

#endif
