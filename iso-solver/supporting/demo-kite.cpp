
#include <complex>
#include <vector>
#include <iostream>
#include <cmath>
using namespace std;

#include "basics.h"


// atih: imaginary part of atan on imaginary axis, atih(y) == Im atan(iy)
double atih(const double x)
{
    if (abs(x)<=1.0) return atanh(x);
    else return atanh(1.0/x);
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
    return 0.5*log( (epsilon*epsilon + delta*delta) );
}



/* calculate the quantum numbers (J) as they are in the complex (non-Takahashi) Bethe equations, from the solution */
vector< complex<double> > dirtyQuantumJforRoots (const int for_chain_length, const vector< complex<double> >& for_roots)
{
    const double threshold_quantum_number = 1e-18;
    
 	vector< complex<double> > bethe_i (for_roots.size());
	for (int alpha=0; alpha < for_roots.size(); ++alpha) {
		/*
		if (abs(real(for_roots[alpha])) < threshold_quantum_number && abs(abs(imag(for_roots[alpha]))-0.5)< threshold_quantum_number ) {
			// numberRoots is number of finite roots. if even, there is no zero; if odd, there is a zero.
			int root_sign = isgn(imag(for_roots[alpha]));
			if (for_roots.size()%2) 
				// odd. there's a zero q.n.
				bethe_i[alpha] = 0.25*(for_chain_length - root_sign*(for_chain_length%4));
			else 
				// even. no zero.
				bethe_i[alpha] = 0.25*(for_chain_length - root_sign*(2 - for_chain_length%4));
			continue;
		}
		*/
		complex<double> lhs = atan(2.0*for_roots[alpha]) * double (for_chain_length);
		for (int beta=0; beta<for_roots.size(); ++beta){ 
			if (alpha!=beta) lhs -= atan(for_roots[alpha] - for_roots[beta]);
		}
		bethe_i[alpha] = lhs/PI;	
	}
	return bethe_i;
}


double betheError (const int for_chain_length, const vector< complex<double> >& for_roots)
{
    vector< complex<double> > qns = dirtyQuantumJforRoots (for_chain_length, for_roots);
    double error = 0.;
    for (int i=0; i<qns.size(); ++i) {
        double rounded = 0.5*round(2.*real(qns[i]));
        error += norm(qns[i] - rounded);
    }
    return error;
}



double eps_step(const double N, const double J, const double J2, const double eps1, const double eps2, const double del0, const double del1) 
{
    // Re eq for real rapidity
    /*
    // this will converge the real rapidity, but the sum of the string will be off. it doesn't converge together with the
    // string equations. with the string-sum equation, the central rapidity wil initially have a nonint quantum number but
    // the sum will be correct. that converges together with the deviation equations. 
    double eps_s =  0.5 * tan( 1./(N-1.) * ( J2*PI   // -0.5*PI 
                        + xi(eps2-eps1, 2.+del1) + xi(eps2-eps1, -del1)
                        + xi(eps1+eps2, 2.+del1) + xi(eps1+eps2, -del1)
                        + xi(eps2, 3.+del0)      + xi(eps2, -1.-del0) 
                        ));
    
    */

    // string center eqn
    /*
        double eps_s =  (1.5 + del1) * tan( 1./N * (  (2.*J+1.+J2)*PI 
                                // scattering with top & bottom roots:
                                + xi(eps1, -del0+del1) + xi(eps1,2.+del0-del1) 
                                + xi(eps2, 3.+del0)    + xi(eps2, -1.-del0)
                                + xi(eps1, 4.+del1+del0) + xi(eps1,-2.-del1-del0)
                                // kinetic terms, excluding the inverted term N*xi(2.*eps1,3.+2.*del1) 
                                - N*xi(2.*eps1,-1.-2.*del1) - N*atan(2.*eps2) 
                                // scattering with opposite string:
                                + xi(2.*eps1, 3.+2.*del1) + xi(2.*eps1, -1.-2.*del1)
                                + 2.*(xi(eps1+eps2, 2.+del1) + xi(eps1+eps, -del1))
                                + atan(2.*eps2) + 2.*atan(2.*eps1)
                                ));
        */
        
        
    return  (0.5) * tan( 1./(N-1.) * (  (2.*J+1.+J2)*PI 
                    // scattering with top & bottom roots:
                    + xi(eps1, -del0+del1) + xi(eps1,2.+del0-del1) 
                    + xi(eps2, 3.+del0)    + xi(eps2, -1.-del0)
                    + xi(eps1, 4.+del1+del0) + xi(eps1,-2.-del1-del0)
                    // kinetic terms with absorbed opposites, excluding the inverted term  
                    - (N-1.)*xi(2.*eps1,-1.-2.*del1) - (N-1.)*xi(2.*eps1,3.+2.*del1) // - (N-1.)*atan(2.*eps2) 
                    // scattering with opposite string:
                    // absorbed + xi(2.*eps1, 3.+2.*del1) + xi(2.*eps1, -1.-2.*del1)
                    + 2.*(xi(eps1+eps2, 2.+del1) + xi(eps1+eps2, -del1))
                    +  2.*atan(2.*eps1) // absorbed: + atan(2.*eps2) 
                    ));
}

double del_step(const double N, const double eps1, const double eps2, const double del0, const double del1) 
{       
        // Im eq:  (N-1) atih (4+2del0) = zeta(eps1, 2+del0+del1) - zeta(eps1, -del0+del1) + ... 
        
        return -2. + 1./(2.*tanh( 1./(N-1.)*(  0.
                                                        - zeta(eps1, -del0+del1)
                            + zeta(eps2, 3.+del0)       - zeta(eps2,-1.-del0)
                            + zeta(eps1, 4.+del0+del1)  //+zeta(eps1, 2.+del0+del1)-zeta(eps1,-2.-del0-del1)
                            )));
}       
 
          
int main() {


//    double eps = 0;
//    double del = 0;
    
    double eps_last = 0.1;
    double del_last = 0.01;
    
    double N=16.0;
   
    
    double eps2 = 0; 
    double del0 = 0;
    
    // these small deviations are extremely important to get the result to converge - wrong sign, no result.
    double eps1 = eps2-1e-10; 
    double del1 = -1e-10;//+1e-20; 
    
    // quantum numbers
    double J = (N-5.)/2.;
    double J2 = -0.5;
    
    
    double eps_s0 = eps_step(N, J, J2, eps1, eps2, del0, del1);
    double del_s0 = del_step(N, eps1, eps2, del0, del1);
        


    for(int i=0; i<60; ++i) {
        
        double eps_s = eps_step(N, J, J2, eps1, eps2, del0, del1);
        double del_s = del_step(N, eps1, eps2, del0, del1);
                            
        // secant method step
        double del_n = ( del0*(del_last-del_s0) - del_last*(del0-del_s) ) /   ( del_last -del_s0 -del0 + del_s );
        //double eps_n = ( eps2*(eps_last-eps_s0) - eps_last*(eps2-eps_s) ) /   ( eps_last -eps_s0 -eps2 + eps_s );

        // normal iteration step
        //del_n = del_s;
        double eps_n = eps_s; 
        
        eps_last = eps2;
        del_last = del0;
        eps_s0 = eps_s;
        del_s0 = del_s;

        // shift both epsilons, as the difference must remain a very small number.
        eps1 += eps_n-eps2;
        eps2 = eps_n;
        del0 = del_n;
                   
        double last_eps1 = eps1;
        double last_del1 = del1;
        // find deviations
        if (true) {
            for(int j=0; j<1;++j) {
                // theta = -xi(eps1-eps2, -del1)
                double theta =  (N-1.)*( xi(2.*eps1,3.+2.*del1) + xi(2.*eps1,-1.-2.*del1) ) - PI*(2.*J+1.) - (
                                   xi(eps1, -del0+del1) + xi(eps1,2.+del0-del1) 
                                 + xi(eps1-eps2, 2.+del1)
                                 + xi(eps1, 4.+del1+del0) + xi(eps1,-2.-del1-del0) 
                                 + xi(eps1+eps2, 2.+del1) + xi(eps1+eps2, -del1) 
                                 + 2.*atan(2.*eps1) 
                                );
                               
                // lnr = -ln(sqrt( sq(eps1-eps2) +sq(del1) ))
                double lnr = (N-1.)*( zeta(2.*eps1,3.+2.*del1) - zeta(2.*eps1,-1.-2.*del1) ) - (
                                   zeta(eps1, -del0+del1) - zeta(eps1,2.+del0-del1) 
                                 + zeta(eps1-eps2, 2.+del1) 
                                 + zeta(eps1, 4.+del1+del0) - zeta(eps1,-2.-del1-del0) 
                                 + zeta(eps1+eps2,2.+del1) - zeta(eps1+eps2, -del1) 
                                 + 2.*atih(2.+2.*del1) 
                              );   
                
                // these are the correct equations, otherwise the qns for the string are off.
                eps1 = eps2 - exp(-lnr)*sin(-theta); 
                del1 = -exp(-lnr)*cos(-theta); 
            }    
        }
        
        double convergence = sq(eps2-eps_last) + sq(del0-del_last);
        double small_convergence = sq(eps1-eps2-last_eps1+eps_last) + sq(del1-last_del1);
        
        cout.precision(4);
        
        cout<<endl;
        cout<<" iter "<<i;
        cout<<" lcnv "<<convergence;
        cout<<" scnv "<<small_convergence;
            
        vector<complex<double> > roots(8);
        roots[0] = I*(2.+del0);
        roots[1] = eps1+I*(1.+del1);
        roots[2] = eps2;
        roots[3] = conj(roots[1]);
        roots[4] = -roots[0];
        roots[5] = -roots[1];
        roots[6] = -roots[2];
        roots[7] = -roots[3];
        
        cout<<" berr "<<betheError(N, roots);
        cout<<endl;    
        
        cout.precision(16);
        cout<<"rts "<<roots<<endl;
        cout.precision(4);
        cout<<"qns  "<<dirtyQuantumJforRoots(N, roots)<<endl;
    }    

    
    return 0;
}




