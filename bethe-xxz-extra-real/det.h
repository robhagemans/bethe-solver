#ifndef DET_H
#define DET_H

#include <vector>
#include "bethe.h"
#include "square.h"

/*************************************************************************************************/
/*************************************************/
/* from Numerical Recipes in C, section 2.1, p. 46:
/* LU decomposition
/* and efficient determinant algorithm.
/*
/**/

extern const char* exc_ZeroArray;
extern const char* exc_Singular;

// LU decomposition
template <class number> void decomposeLU (Square<number>& a, vector<int>& indx, int& d);
template <class number> void backsubLU (Square<number>& a, vector<int>& indx, vector<number>& b);

// Returns the ln of the absolute value of the determinant of matrix a, through LU decomposition
// messes up your matrix in the process. takes less memory and is quicker.
REAL lndetDestroy (Square<REAL>& a);
REAL lndetDestroy (Square< complex<REAL> >& a);

inline REAL lndet (Square<REAL> copy) { return lndetDestroy(copy); };
inline REAL lndet (Square< complex<REAL> > copy) { return lndetDestroy(copy); };


#include "det.def"

#endif
