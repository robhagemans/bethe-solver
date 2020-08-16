/* Rob Hagemans, instituut voor theoretische fysica, universiteit van Amsterdam */

#include <vector>
#include <iostream>
using namespace std;

#include <cmath>

#include <complex>
#define PI 3.141592653589793238462643

const complex<double> I(0,1);

// arctangent: the machine cannot resolve whether a number is on the imaginary axis, but the phase shifts change in that case
inline complex<double> atan(const complex<double> argument)
{
	if (abs(real(argument)) > 1e-18)
		return 0.5*I*( log (1.0 - I*argument) - log (1.0+I*argument) );
	else if (abs(imag(argument)) < 1.0)
		return I*atanh(imag(argument));
	else if (imag(argument) > 0)
		return I*atanh(1.0/imag(argument)) + 0.5*PI;
	else
		return I*atanh(1.0/imag(argument)) - 0.5*PI;
}

template<class number> inline number sgn(const number& x) { return (x<0)?-1:( (x==0)?0:1 ); }
template<class number> inline int isgn(const number& x) { return (x<0)?-1:( (x==0)?0:1 ); }

inline double sq(const double a)
{
	return a*a;
}


inline void swap(double& a, double& b)
{
	const double temp = a;
	a = b;
	b = temp;
}

double modulo (double value, const double period)
{
	while (value < 0) value += period;
	while (value >= period) value -= period;
	return value;
}


const double epsilon = 1e-10;
const double start_distance = 1e-2;
const int max_iter = 1000;


double signmodulo(double value, const double period)
{
//cerr<<value<< "->"<< modulo(value,period) <<":  ";
	double result =  (value = modulo(value,period)) <= period/2.0 ? value : value - period;
	//cerr<<" &"<<result<<"% ";
	return result;
}


/*
inline double indicator2(const double x, const double y)
{
// N =22, @J =19
// N=24, 2J=21
//(-6.17307,2.05713)      (-5.49088,1.1835e-13)   (-5.49088,-1.1835e-13)  (-6.17307,-2.05713) 2J=-47 N=48
	//return atan(2.0*(y+x)) - (0.5*-47.0*PI +  atan(x) + real(atan (y+x - (-6.17307 +2.05713*I)) + atan (y+x - (-6.17307 - 2.05713*I)))   )/48.0;
	return atan(2.0*(y+x)) - (0.5*-47.0*PI +  atan(x) + real(atan (y+x - (-y +2.05713*I)) + atan (y+x - (-y - 2.05713*I)))   )/48.0;
}

inline double indicator1(const double x, const double y)
{
	//double result = atan(2.0*y) - (0.5*-19.0*PI +  atan(-x) )/22.0;
	return atan(2.0*(y)) - (0.5*-47.0*PI +  atan(-x) + real(atan (y - (-y +2.05713*I)) + atan (y - (-y - 2.05713*I) ))   )/48.0;

}
	*/


inline double indicator2(const double x, const double y)
{
	//return atan(2.0*(y+x)) - (0.5*-47.0*PI +  atan(x) + real(atan (y+x - (-6.17307 +2.05713*I)) + atan (y+x - (-6.17307 - 2.05713*I)))   )/48.0;
	return 2.0*(y+x) - tan((0.5*-47.0*PI +  atan(x) + real(atan (y+x - (-6.17307 +2.05713*I)) + atan (y+x - (-6.17307 - 2.05713*I)))   )/48.0);

}

inline double indicator1(const double x, const double y)
{

	//return atan(2.0*(y)) - (0.5*-47.0*PI +  atan(-x) + real(atan (y - (-6.17307 +2.05713*I)) + atan (y - (-6.17307 - 2.05713*I) ))   )/48.0;
	return 2.0*(y) - tan((0.5*-47.0*PI +  atan(-x) + real(atan (y - (-6.17307 +2.05713*I)) + atan (y - (-6.17307 - 2.05713*I) ))   )/48.0);


}




double findRoot1(const double x, const double y_guess = 0.0)
{
	// find the zero for this x

	double y_min = y_guess;
	double y_max = y_guess;
	double f_min, f_max;

	int iter=0;
	// bracket a root
	do {
		y_max += 0.1;
		y_min -= 0.1;

		f_min = indicator1(x, y_min);
		f_max = indicator1(x, y_max);

cout<<y_min<<" "<<f_min<<"   ";
cout<<y_max<<" "<<f_max<<endl;

		if (++iter > max_iter) throw "maxiter";
	} while ( sgn(f_min) == sgn(f_max) );


	// we now have a bracket. zoom in.
	while ( abs(y_max - y_min) >= epsilon ) {
			// interpolate
			//double y_try = (abs(f_min)*y_max + abs(f_max)*y_min) / (abs(f_max)+abs(f_min));

		// bisection is faster for this one
		double y_try = 0.5*(y_min+y_max);
		double f_try = indicator1(x, y_try);
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
cout<<y_min<<" "<<f_min<<"   ";
cout<<y_max<<" "<<f_max<<endl;

		if (++iter > max_iter) throw "maxiter";
	}

	return 0.5*(y_min+y_max);
}




bool findRootForBoth(const double y_trivial_guess, double& x_root, double& y_root)
{

	// first find the trivial root. we only need to check one equation far that.
	double x_trivial = 0.0;
	double y_trivial = findRoot1(x_trivial, y_trivial_guess);

	// now that we have one root, let's look for the other.
	// here we need to check both equations

	// fishing expedition to bracket a second root for positive x
	double x_min = x_trivial + start_distance;
	// find y for this x (stay on root curve for indicator_1), use y_trivial as guess
	double y_min = findRoot1(x_min, y_trivial);
	// get the second indicator function
	double f_min = indicator2(x_min, y_min);

	// outer bracket
	double x_max = x_min;
	double y_max = y_min;
	double f_max;

	int iter = 0;

	// fishing expedition for *one* extra root
	do {
		x_max += 1.0;
		// stay on root curve for indicator_1, use last y as guess
		y_max = findRoot1(x_max, y_max);
		// second indicator function
		f_max = indicator2(x_max, y_max);

cout<<iter<<"  ";
cout<<x_min<<" "<<y_min<<" "<<f_min<<"   ";
cout<<x_max<<"  "<<y_max<<" "<<f_max<<endl;

		if (++iter > max_iter) return false;
	} while ( sgn(f_min) == sgn(f_max) );
cout<<endl;


	double x_try, y_try, f_try;
	// we now have a bracket. zoom in.
	while ( abs(x_max - x_min) + abs(y_max-y_min) >= epsilon ) {
		// interpolate
		//x_try = (abs(f_min)*x_max + abs(f_max)*x_min) / (abs(f_max)+abs(f_min));
		// bisect seems faster again.
		x_try = 0.5*(x_min+x_max);
		// stay on curve. use last y_try as guess
		y_try = findRoot1 (x_try, y_try);
		f_try = indicator2 (x_try, y_try);

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

		if (++iter > max_iter) return false;
	}

cout<<indicator1(x_min, y_min)<<" "<< indicator2(x_min, y_min)<<endl;
cout<<indicator1(x_max, y_max)<<" "<< indicator2(x_max, y_max)<<endl;

	x_root = 0.5*(x_min+x_max);
	y_root = 0.5*(y_min+y_max);
	return true;
}





double findRoot1ForY(const double y, const double x_guess = 0.0)
{
	// find the zero for this y

	double x_min = x_guess;
	double x_max = x_guess;
	double f_min, f_max;

	int iter=0;

	// fishing expedition to bracket a root
	do {
		x_max += 0.1;
		x_min -= 0.1;

		f_min = indicator1(x_min, y);
		f_max = indicator1(x_max, y);
cout<<iter<<"  ";
cout<<x_min<<" "<<f_min<<"   ";
cout<<x_max<<" "<<f_max<<endl;

		if (++iter > max_iter) break;
	} while ( sgn(f_min) == sgn(f_max) );


	// we now have a bracket. zoom in.
	while ( abs(x_max - x_min) >= epsilon ) {
		// interpolate
		double x_try = (abs(f_min)*x_max + abs(f_max)*x_min) / (abs(f_max)+abs(f_min));
		double f_try = indicator1(x_try, y);

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
cout<<iter<<"  ";
cout<<x_min<<" "<<f_min<<"   ";
cout<<x_max<<" "<<f_max<<endl;
	}

	return 0.5*(x_min+x_max);
}



int main()
{
	double x, y;

	// first find the trivial root. we only need to check one equation for that.
	double y_trivial = findRoot1(0.0, -5.49);
cerr<<y_trivial<<endl;
	// then find the root that is aligned with the other rapidities
	double eps = findRoot1ForY(-6.17307, y_trivial+6.17307);
cerr<<eps<<endl;;

	//if (findRootForBoth(2.2, x, y))	cout<< y<<"  "<<y+x<<endl;
}

