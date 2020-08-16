#ifndef SQUARE_H
#define SQUARE_H

#include "exception.h"
#include "basics.h"

#include <string>
using namespace std;

/*************************************************************************************************/
/*************************************************/
/* class Square
/* implementing a square matrix
/*
/* this is JS's SQMat class
/* with some minor changes (naming convention, exceptions)
/*
/**/

extern const string exc_IncompatibleSizes;
extern const string exc_SizeOne;

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


/*************************************************************************************************/
/*************************************************/
/* from Numerical Recipes in C, section 2.1, p. 46:
/* LU decomposition
/* and efficient determinant algorithm.
/*
/**/


//extern const char* exc_ZeroArray;
//extern const char* exc_Singular;


// LU decomposition
template <class number> bool decomposeLU (Square<number>& a, vector<int>& indx, int& d);
template <class number> void backsubLU (Square<number>& a, vector<int>& indx, vector<number>& b);

// Returns the ln of the absolute value of the determinant of matrix a, through LU decomposition
// messes up your matrix in the process. takes less memory and is quicker.
double lndetDestroy (Square<double>& a);
double lndetDestroy (Square< complex<double> >& a);

inline double lndet (Square<double> copy) { return lndetDestroy(copy); };
inline double lndet (Square< complex<double> > copy) { return lndetDestroy(copy); };







template <class number>
Square<number>::Square (int N) : dim(N) , M(new number*[N])
{
	M[0] = new number[N*N];
	for (int i = 1; i < N; i++) M[i] = M[i-1] + N;
}

template <class number>
Square<number>::Square (const Square& rhs) : dim(rhs.dim) , M(new number*[dim])
{
	int i,j;
	M[0] = new number[dim*dim];
	for (i = 1; i < dim; i++) M[i] = M[i-1] + dim;
	for (i = 0; i < dim; i++)
		for (j = 0; j < dim; j++) M[i][j] = rhs[i][j];
}

template <class number>
Square<number>::Square (const number& a, int N) : dim(N) , M(new number*[dim])
{
	int i, j;
	M[0] = new number[dim*dim];
	for (i = 1; i < dim; i++) M[i] = M[i-1] + dim;
	for (i = 0; i < dim; i++) {
		for (j = 0; j < dim; j++) M[i][j] = number(0);  // cast
		M[i][i] = a;
	}
}

template <class number>
Square<number>::Square (const Square& a, const Square& b) : dim (a.dim * b.dim) , M(new number*[a.dim * b.dim])
{
	M[0] = new number[a.dim * b.dim * a.dim * b.dim];
	for (int i = 1; i < a.dim * b.dim; ++i) M[i] = M[i-1] + a.dim * b.dim;

	for (int i1 = 0; i1 < a.dim; ++i1)
	for (int i2 = 0; i2 < a.dim; ++i2)
	for (int j1 = 0; j1 < b.dim; ++j1)
	for (int j2 = 0; j2 < b.dim; ++j2)
		M[i1 * (b.dim) + j1][i2 * (b.dim) + j2] = a[i1][i2] * b[j1][j2];
}


template <class number>
Square<number>::Square (const Square& a, int row_id, int col_id) : dim (a.dim - 1) , M(new number*[dim])
{
	if (dim == 0) throw Exception("Square::Square", exc_SizeOne);

	M[0] = new number[dim * dim];

	for (int i = 1; i < dim; ++i) M[i] = M[i-1] + dim;
	for (int i = 0; i < row_id; ++i) for (int j = 0; j < col_id; ++j) M[i][j] = a[i][j];
	for (int i = row_id; i < dim; ++i) for (int j = 0; j < col_id; ++j) M[i][j] = a[i+1][j];
	for (int i = 0; i < row_id; ++i) for (int j = col_id; j < dim; ++j) M[i][j] = a[i][j+1];
	for (int i = row_id; i < dim; ++i) for (int j = col_id; j < dim; ++j) M[i][j] = a[i+1][j+1];
}



template <class number>
Square<number>& Square<number>::operator= (const Square<number>& rhs)
{
	if (this != &rhs) {
		if (dim != rhs.dim) throw Exception ("Square::operator=", exc_IncompatibleSizes);

		for (int i = 0; i < dim; ++i)
			for (int j = 0; j < dim; ++j)
				M[i][j] = rhs[i][j];
	}
	return *this;
}

template <class number>
Square<number>& Square<number>::operator= (const number& a)
{
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j)
			M[i][j] = number(0);
		M[i][i] = a;
	}
return *this;
}

template <class number>
Square<number>& Square<number>::operator+= (const number& a)
{
	for (int i = 0; i < dim; ++i)
		M[i][i] += a;
	return *this;
}

template <class number>
Square<number>& Square<number>::operator+= (const Square<number>& a)
{
	if (dim != a.dim) throw Exception ("Square::operator+", exc_IncompatibleSizes);
	for (int i = 0; i < dim; ++i)
		for (int j = 0; j < dim; ++j)
			M[i][j] += a[i][j];
	return *this;
}

template <class number>
Square<number>& Square<number>::operator-= (const number& a)
{
	for (int i = 0; i < dim; ++i)
		M[i][i] -= a;
	return *this;
}

template <class number>
Square<number>& Square<number>::operator-= (const Square<number>& a)
{
	if (dim != a.dim) throw Exception ("Square::operator-", exc_IncompatibleSizes);

	for (int i = 0; i < dim; ++i)
		for (int j = 0; j < dim; ++j)
			M[i][j] -= a[i][j];
	return *this;
}

template <class number>
Square<number>& Square<number>::operator*= (const number& a)
{
	for (int i = 0; i < dim; ++i) for (int j = 0; j < dim; ++j) M[i][j] *= a;
	return *this;
}

template <class number>
Square<number>& Square<number>::operator*= (const Square<number>& a)
{
	if (dim != a.dim) throw Exception ("Square::operator*", exc_IncompatibleSizes);

	Square<number> leftarg(*this);  // use copy constructor.
	for (int i = 0; i < dim; ++i)
		for (int j = 0; j < dim; ++j) {
			M[i][j] = 0.0;
			for (int k = 0; k < dim; ++k)
				M[i][j] += leftarg[i][k] * a[k][j];
		}
	return *this;
}


template <class number>
Square<number>::~Square()
{
	if (M != 0) {
		delete[] M[0];
		delete[] M;
	}
}

// outputting a square
template<class number>
ostream& operator<< (ostream& stream, Square<number> out_to_be_put) {

	if (!out_to_be_put.size()) return stream;

	for (int i=0; i< out_to_be_put.size(); ++i)	{
		for (int j=0; j< out_to_be_put.size(); ++j)
			stream << ' ' << out_to_be_put[i][j];
		stream <<endl;
	}

	return stream;
}





/*************************************************************************************************/
/*************************************************/
/* from Numerical Recipes in C, section 2.1, p. 46:
/* LU decomposition
/* and efficient determinant algorithm.
/*
/**/


// LU decomposition
//
template <class number>
	bool decomposeLU (Square<number>& a, vector<int>& indx, int& d)
{
	//const char* here = "decomposeLU";
	const number TINY = 1.0e-20;
	int i, imax, j, k;
	number big, dum, sum, temp;

	int n = a.size();
	if (n == 0) return false; //throw Exception(here, exc_ZeroArray);
	vector< number > vv(n);
	d = 1;
	for (i = 0; i < n; i++) {
		big = 0.0;
		for (j = 0; j < n; j++)
		//      if ((temp = fabs(a[i][j])) > big) big = temp;
			if ((abs(temp = a[i][j])) > abs(big)) big = temp;
		if (big == 0.0) return false; //throw Exception(here, exc_Singular);
		vv[i] = 1.0/big;
	}
	for (j = 0; j < n; j++) {
		for (i = 0; i < j; i++) {
			sum = a[i][j];
			for (k = 0; k < i; k++) sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;
		for (i = j; i < n; i++) {
			sum = a[i][j];
			for (k = 0; k < j; k++) sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			//      if ((dum = vv[i]*fabs(sum)) >= big) {
			if ((abs(dum = vv[i]*sum)) >= abs(big)) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 0; k < n; k++) {
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			d = -d;
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a[j][j] == 0.0) a[j][j] = TINY;
		if (j !=n-1) {
			dum = 1.0/(a[j][j]);
			for (i = j + 1; i < n; i++) a[i][j] *= dum;
		}
	}
	return true;
}

template <class number>
	void backsubLU (Square<number>& a, vector<int>& indx, vector<number>& b)
{
  int i, ii=0, ip, j;
  number sum;

  int n = a.size();
  for (i = 0; i < n; i++) {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if (ii != 0)
      for (j = ii-1; j < i; j++) sum -= a[i][j] * b[j];
    else if (sum != 0.0)
      ii = i + 1;
    b[i] = sum;
  }
  for (i = n - 1; i >= 0; i--) {
    sum = b[i];
    for (j = i + 1; j < n; j++) sum -= a[i][j] * b[j];
    b[i] = sum/a[i][i];
  }
}

#endif
