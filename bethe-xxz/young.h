#ifndef YOUNG_H
#define YOUNG_H

#include <vector>
#include <iostream>
#include "exception.h"

extern const char* exc_Overflow;    // factorial(), choose()
extern const char* exc_InvalidID ;  // young()
extern const char* exc_IndexRange;  // young()

// factorial
long long int factorial (int n);
// binomial coefficient
long long int choose (int upper, int lower);


/** Young tableaux **/

class Young {
public:
	std::vector<int> rows;
	int width;
	int height;
public:
	// set empty, zero-sized tableau
	Young (void) : rows(1,0), width(0), height(1) {};
	// construct young tableau (a vector of its row widths) from id number
	Young (const int max_width, const int max_height, long long int id=0);
	// convert vector<int> to Young
	Young (const int max_width, const int max_height, const std::vector<int>& row_widths);
	// copy
	Young (const Young& original);
	// assign
	Young& operator= (const Young& rhs);
	// destroy
	~Young (void) { };

	// element access: width of a row
	// NOTE: got rid of all overloaded operators here because they seem to confuse the compiler and crash the program.
	// NOTE: this was because it inherited public vector int, I guess...
	int Young::readAt (const int index) const;
	int& Young::writeAt (const int index);
// 	inline int Young::at (const int index) const
// 	{ return read(index); };
// 	inline int Young::operator[] (const int index) const
// 	{ 	return read(index); };
// 	inline int& Young::operator[] (const int index)
// 	{ 	return at(index); }


	// number of tableaux of this width, height
	inline long long int limId(void) const { return choose (width+height, width); };
	// set id and sizes
	void setId(const int max_width, const int max_height, const long long int id);
	// calculate id of young tableau
	long long int id(void) const;
	// transpose this tableau (particle-hole transformation)
	Young& transpose(void);
	// make this a full tableau
	Young& setFull (void);
	// make this an empty tableau
	Young& setEmpty (void);
};

extern const std::vector<Young> NO_SHIFTS;

// outputting a young tableau
std::ostream& operator<< (std::ostream& stream, Young out_to_be_put);


#endif
