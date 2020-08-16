#ifndef BETHE_H
#define BETHE_H

#include <complex>
#include <vector>
#include <iostream>
#include <cmath>

#include <sys/resource.h>
#include <stdio.h>
#include <sstream>
#include <string>
using namespace std;

#include "exception.h"


/** exceptions for function defined here **/
extern const char* exc_ZeroDenominator; //polint
extern const char* exc_NonFiniteArg; // findContinuedFraction

/** some generic definitions **/
typedef double REAL;						// so that we can easily make this e.g. float for speed.
// a big number
#define BOOMBANOOMBA  1e+20

/** machine-dependent stuff **/
// number for which tanh(atanh(x)) == 1.0
#define ATANH_THRESHOLD 18.0

// real part below which the argument of atan is considered to be on the imaginary line, i guess this should be machine_epsilon...
#define EPSILON_ATAN 1e-18

/** some obvious math constants **/
// weird ratio between circumference and diameter of circle
#define PI 3.141592653589793238462643
// squre root of weird ratio between circumference and diameter of circle
extern const REAL SQRT_PI;
// square root of -1
extern const complex<REAL> I;
// infinity. I like compiler warnings.
extern const REAL INFINITE;
// not a number (a rather useless constant, really)
extern const REAL NOT_A_NUMBER;

/** math functions **/

// GCC does define these as macros, but they're not part of the standard:
//template<class number> inline number min(const number &a, const number &b) { return (a)<(b) ? a : b; }
//template<class number> inline number max(const number &a, const number &b) { return (a)>(b) ? a : b; }

// square
template<class number> inline number sq(const number& x) { return x*x; }
// sign
template<class number> inline number sgn(const number& x) { return (x<0)?-1:( (x==0)?0:1 ); }
template<class number> inline int isgn(const number& x) { return (x<0)?-1:( (x==0)?0:1 ); }
// factional part
inline REAL frac(const REAL& argument) { return (argument - trunc(argument)); }
// distance to nearest integer
inline REAL rounderr(const REAL& argument) { return (argument - round(argument)); }

// mod arithmetic: not quite the same as % for negative values
inline int modulo (int value, const int period) { return (value %= period)>=0 ? value : value+period; }

inline REAL modulo (REAL value, const REAL period)
{
	while (value < 0) value += period;
	while (value >= period) value -= period;
	return value;
}


// arctangent
//inline complex<REAL> atan(const complex<REAL> argument)  { return 0.5*I*( log (1.0-I*argument) - log (1.0+I*argument) ) ; }

// arctangent: the machine cannot resolve whether a number is on the imaginary axis, but the phase shifts change in that case
inline complex<REAL> atan(const complex<REAL> argument)
{
	if (abs(real(argument)) > EPSILON_ATAN)
		return 0.5*I*( log (1.0 - I*argument) - log (1.0+I*argument) );
	else if (abs(imag(argument)) < 1.0)
		return I*atanh(imag(argument));
	else if (imag(argument) > 0)
		return I*atanh(1.0/imag(argument)) + 0.5*PI;
	else
		return I*atanh(1.0/imag(argument)) - 0.5*PI;
}

// inverse hyperbolic tangent
//inline complex<REAL> atanh(const complex<REAL> argument) { return 0.5*( log(1.0+argument) - log(1.0-argument) ) ; }


/** tests for non-finite numbers **/

// true if nan
//NOTE: (isnan(x)) and (x != x) had problems, this works.
template<class number> inline bool isNan(number x) { return (!((x<0.0) || (x>=0.0))); }
template<class number> inline bool isNan(complex<number> x) { return (!((real(x)<0.0) || (real(x)>=0.0))) || (!((imag(x)<0.0) || (imag(x)>=0.0))); };

// true if not nan, not inf, not -inf.
template<class number> inline bool finite(number x) { return !isNan(x) && (x != INFINITE) && (x != -INFINITE); }
template<class number> inline bool finite(complex<number> x) { return finite(norm(x)); }


/** stuff that didn't fit elsewhere **/

// continued fraction algorithm
vector<int> findContinuedFraction(REAL REAL_number, REAL precision, int max_length);

// exponential approximation to delta function
inline REAL smoothDelta (const REAL x, const REAL epsilon) { return exp(-sq(x)/(4.0*sq(epsilon)))/(2.0*SQRT_PI*epsilon); }

// polynomial interpolation (fron Numerical Recipes)
void polint (vector<REAL>& xa, vector<REAL>& ya, const REAL x, REAL& y, REAL& dy);

/** UNIX interface and stuff **/

// list directory
vector<string> ls (const string description);
// file exists
bool exists (const string file_name);
// make name unique (by appending .[number] to it)
string uniqueName (const string try_name);
// as in PHP explode
vector<string> explode (const string& str, const string& delimiters = " ");


/** timer **/

// convert seconds to human-readable time string
string humanReadable(double secs);

class Stopwatch {
private:
	rusage *p_usage;	// structure returned by getrusage()
	double start_u; 	// initial user time
	double start_s; 	// initial system time
	double lap_u;		// last lap user time
	double lap_s; 		// last lap system time
	double stop_u;
	double stop_s;
	bool run;
public:
	Stopwatch (bool running=true);
	~Stopwatch(void);

	// reset (and start) stopwatch
	void reset(void);
	// stop
	void stop(void);
	// start
	void start(void);
	// returns time since stopwatch last lap readout in seconds
	double lapUse(void);
	// returns time since stopwatch last lap readout in seconds
	double lapSys(void);
	// returns time since stopwatch reset in seconds
	double getUse(void) const;
	// returns time since stopwatch reset in seconds
	double getSys(void) const;
	// integer number of seconds total time
	inline int seconds(void) const { return (int) round(getUse() + getSys()); };
	// human readable form of time since reset
	inline string humanReadable(void) const { return "user " + ::humanReadable (getUse()) + " sys " + ::humanReadable (getSys()); }
	// human readable form of time since last lap readout
	inline string lapHumanReadable(void) { return "user " + ::humanReadable (lapUse()) + " sys " + ::humanReadable (lapSys()); }

};

/** stream output **/

// separation character for streams
#define SEP '\t'


/** vectors **/

// outputting a vector
template<class number>
ostream& operator<< (ostream& stream, vector<number> out_to_be_put) {
	if (out_to_be_put.size() > 0) stream << out_to_be_put[0];
	if (out_to_be_put.size() > 1)
		for (int i=1; i< out_to_be_put.size(); ++i)
			stream << SEP << out_to_be_put[i];
	return stream;
}



/** nullstream (writes to /dev/null) **/
// thanks to Dietmar Kuehl, see http://www.velocityreviews.com/forums/t279756-null-output-stream.html
struct nullstream: std::ostream {
	nullstream(): std::ios(0), std::ostream(0) {}
 };
// global variable: null stream
extern nullstream CNULL;


#endif
