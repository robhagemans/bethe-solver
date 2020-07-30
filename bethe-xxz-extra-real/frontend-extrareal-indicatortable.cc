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

const double eps = -0.5;

inline double indicator1(const double x, const double y, const double z)
{
	
	double l0 = y;
	double l1 = y+x;
	complex<double> l2 = y+eps+I*z;
	complex<double> l3 = y+eps-I*z;
	return 2.0*(l0) - tan(  (0.5*-47.0*PI +  atan(l0-l1) + real(atan (l0-l2) + atan (l0-l3) ) )/48.0  );
	
}

inline double indicator2(const double x, const double y, const double z)
{
	double l0 = y;
	double l1 = y+x;
	complex<double> l2 = y+eps+I*z;
	complex<double> l3 = y+eps-I*z;
	return 2.0*(l1) - tan(  (0.5*-47.0*PI +  atan(l1-l0) + real(atan (l1-l2) + atan (l1-l3)))/48.0  );

}


inline double indicator3(const double x, const double y, const double z)
{
	
	double l0 = y;
	double l1 = y+x;
	complex<double> l2 = y+eps+I*z;
	complex<double> l3 = y+eps-I*z;
	return real(  2.0*(l2) - tan( (0.5*-45.0*PI +  atan(l2-l0) + atan (l2-l1) + atan (l2-l3))/48.0) );
	
}


inline double indicator4(const double x, const double y, const double z)
{
	
	double l0 = y;
	double l1 = y+x;
	complex<double> l2 = y+eps+I*z;
	complex<double> l3 = y+eps-I*z;
	return imag(  2.0*(l2) - tan( (0.5*-45.0*PI +  atan(l2-l0) + atan (l2-l1) + atan (l2-l3))/48.0) );
	
}

int main()
{
	double x, y, z;
		
	double xdiv = 1280.0;
	double ydiv = 1280.0;
		
	cout.precision(5);
	//cout<<"Table,";
	
	
	for (int iz= -5; iz<5; ++iz) {
		z = 7.5375 + iz /160.0;
		
		cout<<endl<<endl;
		cout<<z<<endl;
		//for (int iy =-80; iy<0; ++iy) {
		//	y = iy/ydiv;
		//	if (y>0) cout<<"+"; else if (y<0) cout<<"-"; else cout<<"0";
		//}
		//cout<<endl;	
		cout<<endl;
	
	
		for (int ix =-20; ix<20; ++ix) {
			x = 2.836 + ix/xdiv;
			if (x>0) cout<<"+ "; else if (x<0) cout<<"- "; else cout<<"0 ";
			//cout<<x<<",";
			for (int iy =-40; iy<40; iy++) {
				y = -6.919 + iy/ydiv;
				//cout<<sq(indicator2(x,y)) + sq(indicator1(x,y))<<",";
				char c;
				if (indicator2(x,y,z) >0 && indicator1(x,y,z)>0) c = '.';
				if (indicator2(x,y,z) <0 && indicator1(x,y,z)<0) c = '\"';
				if (indicator2(x,y,z) <0 && indicator1(x,y,z)>0) c = 'x';
				if (indicator2(x,y,z) >0 && indicator1(x,y,z)<0) c = '=';
				
				cout<<c;
				
				if (indicator3(x,y,z)>0) c = '+'; else if (indicator3(x,y,z)<0) c='-'; else c='0';
				cout<<c;
			}
			cout<<endl;
		}
	}		
	
	x= 2.836;
	y=-6.919;
	z=7.538;
	//eps=-0.5;
	cout<<indicator1(x,y,z)<<" "<<indicator2(x,y,z)<<" "<<indicator3(x,y,z)<<" "<<indicator4(x,y,z)<<endl;
}



