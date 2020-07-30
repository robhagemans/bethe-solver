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


complex<double> q(const complex<double> l0, const complex<double> l1, const complex<double> l2, const complex<double> l3)
{
		
	return (48.0*atan(2.0*l0) - atan(l0-l1) - atan (l0-l2) - atan (l0-l3) )/PI ;

}



int main()
{
    complex<double> l0 (-6.17307, 2.05713);
    complex<double> l1 (-6.17307, 2.05713);
    complex<double> l2 (-4.547135, 0.0);
    complex<double> l3 (-2.8036 , 0.0); //-6.2906
    
    cout << q(l0,l1,l2,l3) <<endl;
    cout << q(l1,l0,l2,l3) <<endl;
    cout << q(l2,l1,l0,l3) <<endl;
    cout << q(l3,l1,l2,l0) <<endl;

    return 0;
}

