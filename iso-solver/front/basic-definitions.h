#ifndef BASIC_DEFS_H
#define BASIC_DEFS_H

#include <complex>
#include <vector>
#include <iostream>
#include <cmath>


/* math constants */

const double PI = M_PI;
const std::complex<double> I (0.0, 1.0);		



/* the square */

// square
template<class number> 
inline number sq(const number& x) 
{ 	return x*x; }



/* signs */


// integer sign
template<class number> 
inline int isgn(const number& x) 
{ 	return (x<0)?-1:( (x==0)?0:1 ); }

// integer sign with sign(0) = +1
template<class number> 
inline int isgn_plus(const number& x) 
{ 	return (x<0)? (-1): (1); }

// integer sign with sign(0) = -1
template<class number> 
inline int isgn_minus(const number& x) 
{ 	return (x<=0)? (-1): (1); }


template<class number> 
inline number sgn(const number& x) 
{	return (number) isgn(x); }

template<class number> 
inline number sgn_plus(const number& x) 
{	return (number) isgn_plus(x); }


/* vectors */

// outputting a vector
template<class number>
std::ostream& operator<< (std::ostream& stream, std::vector<number> out_to_be_put) {
	stream<<"[";
	if (out_to_be_put.size() > 0) { 
		for (int i=0; i< out_to_be_put.size(); ++i) {
	        if (i>0) stream<<' ';		
			stream << out_to_be_put[i];
		}	
    }
	stream<<"]";
	
	return stream;
}



/* empty functional */
// used in SimpleScanner, Solver
class NoFunc {
public:
    inline void operator()(...) {};
};

#endif
