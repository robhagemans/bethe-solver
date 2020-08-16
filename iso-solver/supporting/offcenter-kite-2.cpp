#include <iostream>
#include "basics.h"

using namespace std;





// atih: imaginary part of atan on imaginary axis, atih(y) == Im atan(iy)
double atih(const double x)
{
    if (abs(x)<=1.0) return atanh(x);
    else return atanh(1.0/x);
}

// Re atan iy
double rati(const double y)
{
    if (abs(y)<1.0) return 0.;
    if (y==1.0) return 0.25*PI; //0.25
    if (y==-1.0) return -0.25*PI; //0.25
    return 0.5*PI*sgn(y);
}


/** complex atan sums and differences **/

// argument-like function or 'continuous arctan', see notes.
// atan(a+ib) + atan(a-ib) = xi(a, 1+b) + xi(a,1-b)
inline double xi (const double epsilon, const double delta) 
{
    return (atan(epsilon/delta) + ((delta<0.0)?PI*sgn(epsilon):0));
}	


// atan(a+ib) - atan(a-ib) = i (zeta(a,1+b) - zeta(a,1-b))
inline double zeta(const double epsilon, const double delta) 
{
    return log( sqrt(epsilon*epsilon + delta*delta) );
}


template <typename number>
number secantStep(
        const number z_last, const number z, 
        const number s_last, const number s, 
        double& convergence)
{ 
        number z_next = ( z*(z_last-s_last) - z_last*(z-s) ) /   ( z_last -s_last -z +s );
        // secant method convergence
        convergence += norm(z_next - z);
        return z_next;
}

template <typename number>
number dampingStep(
        const number z, 
        const number s, 
        const double damping,
        double& convergence)
{ 
        number z_next = damping*s + (1.0-damping)*z;
        convergence += norm(z_next - z);
        return z_next;
}

/** Kite class **/


     
    
    
vector< complex<double> > step(const double big_n, const vector<complex<double> >& roots) 
{
    vector< complex<double> > scatter(4, 0.);
    vector< complex<double> > new_roots(4, 0.);
    
    for (int i=0; i<3;++i) {
        for (int j=0; j<4;++j) {
            if (j==i) continue;
            scatter[i] += atan(roots[i]-roots[j]);
        }
    }
    scatter[0] += -0.5*PI;
    scatter[1] += -1.5*PI;
    scatter[2] += 1.5*PI;
    
    
    
    double big_lam1 = tan(real(scatter[1])/big_n ); 
    double eps1 = 2.*big_lam1/(1.-2.*sqrt(0.75)*big_lam1);
    
    double big_lam2 = tan(real(scatter[2])/big_n  ); 
    double eps2 = 2.*big_lam2/(1.-2.*sqrt(0.75)*big_lam2);
    
    
    double del0 = imag(roots[0])-1.;
    double big_del0 = real(roots[0]) - sqrt(0.75);
    
    complex<double> idel0mbigdel0 = tan( scatter[0]/big_n - 0.25*I*log(3.) ) 
                        * 2.*( sqrt(0.75)*big_del0 + del0 + I*(big_del0 + (2.+del0)*sqrt(0.75)) );
         
    new_roots[0] = sqrt(0.75) + I +idel0mbigdel0; //-conj(
    new_roots[1] = sqrt(0.75) + eps1;
    new_roots[2] = sqrt(0.75) + eps2;     
    new_roots[3] = conj(new_roots[0]);
    
    return new_roots;
}

vector< complex<double> > z_last_;
vector< complex<double> > s_last_;



bool iterate(const double big_n, vector< complex<double> >& roots, double& convergence)
{
    
    if (!z_last_.size()) {
        z_last_ = roots;
        s_last_ = step(big_n, roots);
        roots = s_last_;
    }
    vector< complex<double> >  s = step(big_n, roots);
    vector< complex<double> >  z_next(4);
    for (int i=0; i<4; ++i) {
        z_next[i] = secantStep(z_last_[i], roots[i], s_last_[i], s[i], convergence);
    }
    z_next[1] = s[1];
    z_next[2] = s[2];
    z_last_ = roots;
    s_last_ = s;

    for (int i=0; i<4;++i) {
        convergence  += norm(roots[i] - s[i]);
    
    }
    roots = z_next;
    
    return true;
}



/* calculate the quantum numbers (J) as they are in the complex (non-Takahashi) Bethe equations, from the solution */
vector< complex<double> > dirtyJ (const int for_chain_length, const vector< complex<double> >& for_roots)
{
    //const double threshold_quantum_number = 1e-18;
    
 	vector< complex<double> > bethe_i (for_roots.size());
	for (int alpha=0; alpha < for_roots.size(); ++alpha) {
		complex<double> lhs = atan(2.0*for_roots[alpha]) * double (for_chain_length);
		for (int beta=0; beta<for_roots.size(); ++beta){ 
			if (alpha!=beta) lhs -= atan(for_roots[alpha] - for_roots[beta]);
		}
		bethe_i[alpha] = lhs/PI;	
	}
	return bethe_i;
}


int main()
{
    vector< complex<double> > roots = { sqrt(0.75)+1.5*I, sqrt(0.75)+1. , sqrt(0.75)-1., sqrt(0.75)-1.5*I };
//cerr<<roots.size()<<endl;
    for (int iter=0; iter<1000; ++iter) {
 
        double conv = 0.;
        iterate(12, roots, conv);
 
        cout<<iter<<" "<<conv<<":  "<<roots<<endl;
 
    }    
    cout<<" :  "<<roots<<endl;
 
    cout<<dirtyJ(12, roots)<<endl;
}
