#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include "exception.h" 

extern const char* exc_IncompatibleHeights; //Matrix class
extern const char* exc_IncompatibleWidths; //Matrix class
extern const char* exc_NoStream; //Matrix class
extern const char* exc_IndexRange; 

/** non-square matrix (just a 2D lot of numbers really) **/

/// TODO: this should be implemented somewhat more seriously
template<class number>
class Matrix: public vector< vector< number> > {
public:
	Matrix (int height, int width);
	// adding matrices
	Matrix<number>& operator+= (const Matrix<number>& argument);
	// scalar multiplication
	Matrix<number>& operator*= (const number argument);
	// width
	inline int width() const {  return vector< vector< number> >::at(0).size();  };
	// height
	inline int height() const {  return vector< vector<number> >::size();  };
};

// outputting a matrix
template<class number> ostream& operator<< (ostream& stream, Matrix<number> out_to_be_put);





/** template implementations **/

/** class Matrix **/
/** this is a quick-and-dirty implementation. as soon as we really use this, streamline this a bit **/

template<class number>
Matrix<number>::Matrix (int height, int width) 
	: vector< vector<number> > (height, vector<number>(width))
{
}

template<class number>
Matrix<number>& Matrix<number>::operator+= (const Matrix<number>& argument)
{
	const char* here = "Matrix::operator+=";
	
	if (argument.height() != height()) throw Exception (here, exc_IncompatibleHeights);
	if (!height()) return *this; // both empty matrices, do nothing, no problem.
	if (argument.width() != width()) throw Exception (here, exc_IncompatibleWidths);
	for (int i=0; i< height(); ++i)
	for (int j=0; j< width(); ++j) 
		(*this)[i][j] += argument[i][j];
	return *this;
}


template<class number>
Matrix<number>& Matrix<number>::operator*= (const number argument)
{
	for (int i=0; i< height(); ++i)
	for (int j=0; j< width(); ++j) 
		(*this)[i][j] *= argument;
	return *this;
}


// outputting a matrix
template<class number>
ostream& operator<< (ostream& stream, Matrix<number> out_to_be_put) {
	
	if (!out_to_be_put.size()) return stream;
	
	for (int i=0; i< out_to_be_put.height(); ++i)	{
		for (int j=0; j< out_to_be_put.width(); ++j)	
			stream << SEP << out_to_be_put[i][j];
		stream <<endl;
	}
	
	return stream;
}

template<class number>
istream& operator>> (istream& stream, Matrix<number> in_to_be_put)
{
	const char* here= "operator>> into Matrix";
	if (!in_to_be_put.size()) return stream;
	
	for (int i=0; i< in_to_be_put.size(); ++i)	{
		for (int j=1; j< in_to_be_put[0].size(); ++j)	{
			if (!stream) throw Exception (here, exc_NoStream);
			stream >> in_to_be_put[i][j];
		}
	}
	return stream;
}


#endif
