/* Rob Hagemans, instituut voor theoretische fysica, universiteit van Amsterdam */

#include <vector>
#include <iostream>
using namespace std;

#include <cmath>
#define PI 3.141592653589793238462643

template<class number> inline number sgn(const number& x) { return (x<0)?-1:( (x==0)?0:1 ); }
template<class number> inline int isgn(const number& x) { return (x<0)?-1:( (x==0)?0:1 ); }



inline void swap(double& a, double& b)
{
	const double temp = a;
	a = b;
	b = temp;
}



const double epsilon = 1e-10;
const int max_iter = 100;


inline double indicator_2(const double x, const double y)
{
	return 22.0 * atan(2.0*y+x) - ( 0.5*PI* 19.0 +  atan(x) );
}

inline double indicator_1(const double x, const double y)
{
	return 22.0 * atan(2.0*y) - ( 0.5*PI* 19.0 +  atan(-x) );
}



double find_zero_1(const double x, const double y_guess = 0.0)
{
	// find the zero for this x

	double y_min = y_guess;
	double y_max = y_guess;
	double f_min, f_max;

	int iter=0;
	// bracket a root
	do {
		y_max += 1.0;
		y_min -= 1.0;

		f_min = indicator_1(x, y_min);
		f_max = indicator_1(x, y_max);

		if (++iter > max_iter) throw "maxiter";
	} while ( sgn(f_min) == sgn(f_max) );


	// we now have a bracket. zoom in.
	while ( abs(y_max - y_min) >= epsilon ) {
			// interpolate
			//double y_try = (abs(f_min)*y_max + abs(f_max)*y_min) / (abs(f_max)+abs(f_min));

		// bisection is faster for this one
		double y_try = 0.5*(y_min+y_max);
		double f_try = indicator_1(x, y_try);
		if (sgn(f_try)==0) {
			// not bloody likely
			y_min = y_max = y_try;
			f_min = f_max = f_try;
		}
		else if (sgn(f_try) == sgn(f_min)) {
			y_min = y_try;
			f_min = f_try;
		}
		else {
			y_max = y_try;
			f_max = f_try;
		}
		if (++iter > max_iter) throw "maxiter";
	}

	return 0.5*(y_min+y_max);
}




int main()
{
	//double y = 2.2;

	// first find the trivial root. we only need to check one equation far that.
	double x_trivial = 0;
	double y_trivial = find_zero_1(x_trivial, 2.2);

	// now that we have one root, let's look for the other.
	// here we need to check both equations

	// fishing expedition to bracket a second root
	// we have inner and outer brackets
	double x_min_in = x_trivial - epsilon;
	double x_max_in = x_trivial + epsilon;
	// find y for this x (stay on root curve for indicator_1), use y_trivial as guess
	double y_min_in = find_zero_1(x_min_in, y_trivial);
	double y_max_in = find_zero_1(x_max_in, y_trivial);
	// get the second indicator function
	double f_min_in = indicator_2(x_min_in, y_min_in);
	double f_max_in = indicator_2(x_max_in, y_max_in);

	// outer brackets
	double x_min_out = x_min_in;
	double x_max_out = x_max_in;
	double y_min_out = y_min_in;
	double y_max_out = y_max_in;

	int iter = 0;
	double f_min_out, f_max_out;


	// fishing expedition for *one* extra root
	// NOTE: we may want to change y instead of x -- might converge faster. check this.
	do {
		x_max_out += 1.0;
		x_min_out -= 1.0;

		// stay on root curve for indicator_1, use last y as guess
		y_max_out = find_zero_1(x_max_out, y_max_out);
		y_min_out = find_zero_1(x_min_out, y_min_out);

		// second indicator function
		f_min_out = indicator_2(x_min_out, y_min_out);
		f_max_out = indicator_2(x_max_out, y_max_out);

cout<<iter<<"  ";
cout<<x_min_out<<" "<<y_min_out<<" "<<f_min_out<<"   ";
cout<<x_min_in<<"  "<<y_min_in<<" "<<f_min_in<<"   ";
cout<<x_max_in<<"  "<<y_max_in<<" "<<f_max_in<<"   ";
cout<<x_max_out<<"  "<<y_max_out<<" "<<f_max_out<<endl;

		if (++iter > max_iter) throw "maxiter";
	} while ( (sgn(f_min_out) == sgn(f_min_in)) && (sgn(f_max_out) == sgn(f_max_in)) );
cout<<endl;


	// which side is our extra root on?
	double x_min, x_max, y_min, y_max, f_min, f_max;
	/*
	if (sgn(f_max_out) == sgn(f_max_in))
	 {
		// we're bracketing on the min side
		x_min = x_min_out;
		x_max = x_min_in;

		y_min = y_min_out;
		y_max = y_min_in;

		f_min = f_min_out;
		f_max = f_min_in;
	} else
	*/
	// we always take positive side.
	// TODO: get rid of fishing on the negative side if we're going to use x as parameter value
	{
		// we're bracketing on the max side
		x_min = x_max_in;
		x_max = x_max_out;

		y_min = y_max_in;
		y_max = y_max_out;

		f_min = f_max_in;
		f_max = f_max_out;
	}


	double y_try;
	// we now have a bracket. zoom in.
	while ( abs(x_max - x_min) + abs(y_max-y_min) >= epsilon ) {
		// interpolate
		double x_try = (abs(f_min)*x_max + abs(f_max)*x_min) / (abs(f_max)+abs(f_min));
		// stay on curve. use last y_try as guess
		y_try = find_zero_1 (x_try, y_try);
		double f_try = indicator_2(x_try, y_try);

		if (sgn(f_try)==0) {
			// not bloody likely
			x_min = x_max = x_try;
			y_min = y_max = y_try;
			f_min = f_max = f_try;
		}
		else if (sgn(f_try) == sgn(f_min)) {
			x_min = x_try;
			y_min = y_try;
			f_min = f_try;
		}
		else {
			x_max = x_try;
			y_max = y_try;
			f_max = f_try;
		}
cout<<iter<<"  ";
cout<<x_min<<" "<<y_min<<" "<<f_min<<"   ";
cout<<x_max<<"  "<<y_max<<" "<<f_max<<endl;

		if (++iter > max_iter) throw "maxiter";
	}


cout<<indicator_1(x_min, y_min)<<" "<< indicator_2(x_min, y_min)<<endl;
cout<<indicator_1(x_max, y_max)<<" "<< indicator_2(x_max, y_max)<<endl;


	return 0;
}



/*
double find_zero_1_for_y(const double y)
{
	// find the zero for this y

	double x_min = 0.0;
	double x_max = 0.0;
	double f_min, f_max;

	int iter=0;

	// fishing expedition to bracket a root
	do {
		x_max += 1.0;
		x_min -= 1.0;

		f_min = indicator_1(x_min, y);
		f_max = indicator_1(x_max, y);

		if (++iter > max_iter) break;
	} while ( sgn(f_min) == sgn(f_max) );


	// we now have a bracket. zoom in.
	while ( abs(x_max - x_min) >= epsilon ) {
		// interpolate
		double x_try = (abs(f_min)*x_max + abs(f_max)*x_min) / (abs(f_max)+abs(f_min));
		double f_try = indicator_1(x_try, y);

		if (sgn(f_try)==0) {
			// not bloody likely
			x_min = x_max = x_try;
			f_min = f_max = f_try;
		}
		else if (sgn(f_try) == sgn(f_min)) {
			x_min = x_try;
			f_min = f_try;
		}
		else {
			x_max = x_try;
			f_max = f_try;
		}
		if (++iter > max_iter) break;
	}

	return 0.5*(x_min+x_max);
}
*/
