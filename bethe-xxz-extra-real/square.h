#ifndef SQUARE_H
#define SQUARE_H

#include "exception.h"
#include "bethe.h"

/*************************************************************************************************/
/*************************************************/
/* class Square
/* implementing a square matrix
/*
/* this is JS's SQMat class
/* with some minor changes (naming convention, exceptions)
/*
/**/


extern const char* exc_IncompatibleSizes;
extern const char* exc_BaldChicken;

template <class number>
class Square {
private:
	int dim;
	number** M;
public:
	Square (int N);             			// initializes all elements of this n by n matrix to zero
	Square (const Square& rhs);  			// copy constructor
	Square (const number& a, int N);  			// initialize to diagonal matrix with value a (NOT like in NR !!!)
	Square (const Square& a, const Square& b); 	// initialize to tensor product of a and b
	Square (const Square& a, int row_id, int col_id); // init by cutting row row_id and col col_id
	~Square();

	inline number* operator[] (const int i) { return M[i]; };  		// subscripting: pointer to row i
	inline const number* operator[] (const int i) const { return M[i]; };
	inline int size() const { return dim; };

	Square& operator= (const Square& rhs);  		// assignment
	Square& operator= (const number& a);        		// assign 1 to diagonal elements (NOT like in NR !!!)
	Square& operator+= (const number& a);
	Square& operator+= (const Square& a);
	Square& operator-= (const number& a);
	Square& operator-= (const Square& a);
	Square& operator*= (const number& a);
	Square& operator*= (const Square& a);

};

#include "square.def"
#endif
