
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
    return log( sqrt(epsilon*epsilon + delta*delta) );
}


/* string roots */
struct roots {
    vector<double> epsilon; // real shift on either side of the imaginary axis
    vector<double> delta; // imaginary shift from string position    
};

/*

vector< complex<double> > get_roots(const roots& in)
{
    
    vector< complex<double> > out;
    
    return out;
}

*/

/* calculate the quantum numbers (J) as they are in the complex (non-Takahashi) Bethe equations, from the solution */
vector< complex<double> > dirtyQuantumJforRoots (const int for_chain_length, const vector< complex<double> >& for_roots)
{
    const double threshold_quantum_number = 1e-18;
    
 	vector< complex<double> > bethe_i (for_roots.size());
	for (int alpha=0; alpha < for_roots.size(); ++alpha) {
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
		complex<double> lhs = atan(2.0*for_roots[alpha]) * double(for_chain_length);
		for (int beta=0; beta<for_roots.size(); ++beta){ 
			if (alpha!=beta) lhs -= atan(for_roots[alpha] - for_roots[beta]);
		}
		bethe_i[alpha] = lhs/PI;	
	}
	return bethe_i;
}



/** calculate a new value for the deviance (delta) and aberration (epsilon) **/
roots step (const int chain_length, const roots& roots_in)
{
	roots roots_out = roots_in;

cout<< roots_out.epsilon  <<endl;
cout<< roots_out.delta <<endl;
	
	// these are for odd string lengths:
	double big_n_m1 = chain_length -1;
    int length_j = roots_in.epsilon.size()*2 - 1;
    int number_roots = roots_in.epsilon.size()*4 - 4;
    
    // 'small letter' thetas and l's. 
    // these are the rhs of the conjugate-sum and conjugate-diff equation.
	// in the notes, l_0 = ln r_0 = Im(equation_0)    
	vector<long double> theta(length_j/2+1);
	vector<long double> log_rsq(length_j/2+1);
	
	for (int a=0; a<=length_j/2; ++a) {
		double del = roots_in.delta[a]; 
		double eps = roots_in.epsilon[a];

	    theta[a] = 0.0;
    	log_rsq[a] = 0.0;
	
	
	    // ** quantum number contribution to theta **
	    // we don't actually need to know the quantum numbers.
		// contribution is M mod 2 from sum of conjugate quantum numbers
		// this just means that J_+ + J_- = 2J+1  (except for an inner two-pair, which may be a narrow pair)
		// we're adding \pi(2J+1) but since we're taking a cos or sin at the end, all that matters is whether it's even or odd.
		// Js are int for M odd, so 2J+1 odd
		// Js are half-int for M even, so 2J+1 even
		//theta -= PI*number_roots;		
		if (a!=length_j/2) {
		    if (number_roots%2) theta[a] -= PI;
		}
		else {
		    // real equation - eq+ eq^* = 2 Re(eq)
		    // so qn contrib is 2J. assume J=1/2 for M even, J=1 for M odd:
		    // we need to know mod 4 here, as we're doing sin or cos of theta/2.
		    theta[a] += PI;
		    if (number_roots%2) theta[a] += PI;
		}
		
		// ** kinetic term **
		// includes scattering with the opposite root through (N-1)
		
		// lambda = eps + i ( (l_j-1)/2 - a + delta)  
		
		theta[a] += big_n_m1 * ( 
				  xi(2.0*(eps), length_j - 2.0*a + 2.0*del) 
				+ xi(2.0*(eps), -length_j + 2.0*(a+1) - 2.0*del)
			);
		log_rsq[a] += big_n_m1 * (
				   log(sq(eps) + sq(0.5*(length_j) - (a+1) + del))
				 - log(sq(eps) + sq(0.5*(length_j) - (a) + del))
			);
cout<<"1"<<endl;
cout<< theta  <<endl;
cout<< log_rsq <<endl;


        // ** scattering with self-conjugate roots, if applicable **
        // for odd strings, self-conjugate terms are non-dangerous.
        // for odd strings, root 0 and root n/2 are pure real and pure imag, hence excluded here.
        if (a!=0 && a!=length_j/2) {
            // scattering with conjugate opposite, a pure real term, counted twice as we're doing eq_0 + eq_0*
            theta[a] -= 2.0* atan (2.0*eps);
            // scattering with conjugate - a pure imaginary term, also counted twice as we're doing eq_0 - eq_0*
            // in addition, our log_rsq counts all the other terms twice, too, so another factor 2 here.
            log_rsq[a] -= 4.0* atih (length_j-1 -2*a + 2.0*del);
        } 
cout<<"2"<<endl;        
cout<< theta  <<endl;
cout<< log_rsq <<endl;


        // ** scattering with pure imaginary top and bottom roots **
        if (a!=0) {
			double del_0 = roots_in.delta[0]; 
            if (a!=1) {
                // exclude terms if dangerous
                theta[a]    -= xi( eps, 1.0-a+(del-del_0) );
			    log_rsq[a]  += log( sq(eps) + sq(1.0-a+(del-del_0)) );
		    }
		    if (a!=length_j-2) {
		        theta[a]    -= xi( eps, 2.0-length_j+a-(del+del_0) );
                log_rsq[a]  -= log( sq(eps) + sq(2.0-length_j + a - (del+del_0)) );
            }
		    theta[a] -= xi( eps, 1.0+a-(del-del_0) );
            theta[a] -= xi( eps, length_j-a+(del+del_0) );
            
		    log_rsq[a] -=  log (sq(eps) + sq(1.0+a-(del-del_0)) );
            log_rsq[a] +=  log (sq(eps) + sq(length_j-a+(del+del_0)) );
        }
        
cout<<"3"<<endl;        
cout<< theta  <<endl;
cout<< log_rsq <<endl;

		// ** scattering with the rest of the string complex, excluding dangerous terms **		
		for (int b=1; b < length_j-1; ++b) {
			if (b==a) continue;
			double del_b = roots_in.delta[b]; 
			double eps_b = roots_in.epsilon[b];
			
			if (a-b != 1) {
			    theta[a] -= xi( (eps-eps_b), (1.0-(a-b)+(del-del_b)) );
			    theta[a] -= xi( (eps+eps_b), (1.0-(a-b)+(del-del_b)) );
			
			    log_rsq[a] += log( sq(eps-eps_b)  +  sq(1.0-(a-b)+(del-del_b)) );
			    log_rsq[a] += log( sq(eps+eps_b)  +  sq(1.0-(a-b)+(del-del_b)) );
			}
			if (b-a != 1) {
			    theta[a] -= xi( (eps-eps_b), (1.0+(a-b)-(del-del_b)) );
   			    theta[a] -= xi( (eps+eps_b), (1.0+(a-b)-(del-del_b)) );

    			log_rsq[a] -= log( sq(eps-eps_b)  +  sq(1.0+(a-b)-(del-del_b)) );
                log_rsq[a] -= log( sq(eps+eps_b)  +  sq(1.0+(a-b)-(del-del_b)) );
            }
		}		
 
cout<<"4"<<endl;        
cout<< theta  <<endl;
cout<< log_rsq <<endl;

		// ** scattering with other roots **
		// not implemented yet
	}		
cout<< theta  <<endl;
cout<< log_rsq <<endl;

    // sum equations 
    // these equal normsqs and arguments of product vectors
    vector<long double> big_theta(length_j/2+1);
	vector<long double> big_log_rsq(length_j/2+1);
	big_log_rsq[0] = log_rsq[0];  //?? 0.5* // above, we've calculated Im(eq_0 - eq_0^*) rather than Im(eq_0) as in notes. 
	big_theta[length_j/2] = theta[length_j/2];
	for (int a=1; a<=length_j/2; ++a) {
		big_log_rsq[a] = big_log_rsq[a-1] + log_rsq[a];
		big_theta[length_j/2-a] = big_theta[length_j/2-a+1] + theta[length_j/2-a];
    }	
	// we should have big_theta[0]=0 ?
cout<< big_theta  <<endl;
cout<< big_log_rsq <<endl;

cout<<" --> "<<cos(big_theta[length_j/2-1])<< " " << sin(big_theta[length_j/2-1]) <<endl;


cout<< "eps: "<< exp(0.25*big_log_rsq[0])*sin(PI-0.5*big_theta[1]) <<endl;
cout<< "del: "<< exp(0.25*big_log_rsq[0])*cos(PI-0.5*big_theta[1]) <<endl;

    // derive epsilon and delta from product vectors
    roots_out.epsilon[0] = 0;
    roots_out.delta[length_j/2] = 0;
    
    
    roots_out.delta[length_j/2-1] = sqrt( 0.5* ( 
        exp(0.5*big_log_rsq[length_j/2-1])*cos(big_theta[length_j/2-1]) + exp(0.5*big_log_rsq[length_j/2])*cos(big_theta[length_j/2]) 
        ));
    roots_out.epsilon[length_j/2] = 0.5*exp(0.5*big_log_rsq[length_j/2])*sin(big_theta[length_j/2]) / roots_out.delta[length_j/2-1];
 
 /*       
    for (int a=length_j/2; a>0; --a) {
        long double r_am1 = exp(0.5*big_log_rsq[a-1]);
        long double rcosth_am1 = r_am1 * cos(big_theta[a-1]);
        double eps_a = roots_out.epsilon[a];
        if (a>1) 
            roots_out.epsilon[a-1] = sqrt(  0.5*( eps_a - rcosth_am1 + sqrt(sq(sq(eps_a)) - eps_a*rcosth_am1 + sq(r_am1)) )  );
        if (a<length_j/2) 
            roots_out.delta[a-1] = roots_out.delta[a] + 0.5*r_am1*sin(big_theta[a-1])/ roots_out.epsilon[a-1]; 		
	}
	*/
	
		
cout<< roots_out.epsilon  <<endl;
cout<< roots_out.delta <<endl;
exit(0);

    return roots_out;
}





vector<double> secant_step(
        const vector<double>& z_last, const vector<double>& z, 
        const vector<double>& s_last, const vector<double>& s, 
        double& convergence)
{ 
        vector<double> z_next(z.size());
        for (int i = 0; i < z_next.size(); ++i) {
            // secant method step
            z_next[i] = ( z[i]*(z_last[i]-s_last[i]) - z_last[i]*(z[i]-s[i]) ) /   ( z_last[i] -s_last[i] -z[i] +s[i] );
            // secant method convergence
            convergence += (z_next[i] - z[i])*(z_next[i] - z[i]);
        }
        return z_next;
}

vector<double> simple_step(
        const vector<double>& z, 
        const vector<double>& s, 
        const double damping,
        double& convergence)
{ 
        vector<double> z_next (z.size());
        for (int i = 0; i < z_next.size(); ++i) {
            z_next[i] = damping*s[i] + (1.0-damping)*z[i];
            convergence += (z_next[i] - z[i])*(z_next[i] - z[i]);
        }
        return z_next;
}



int main()
{
    // N atih 2y == sum_y+ atih (y - yn) + atih (y + yn)
    // atih inverse is two valued, tanh(x) for narrow strings, 1/tanh(x) for wide strings
    cout.precision(15);


    // single kite
    //#3#0: 2I [ +0 ;; +0 ;]  2J [-1 +9 +1 -9] sums [-1;+1;] sum +0 (+0.01414495,+0) (+0,+1.0035739) (-0.01414495,+0) (+0,-1.0035739)
    //iter 8  (  ), i(  )  ( 0.0141449501592471 ), i( 1.00357387588898 )  conv 2.27585135999056e-21
    int N = 10;

    roots z_last = {  {0.0, 0.014},  {0.0034, 0.0} };
    roots z = {  {0.0, 0.0141}, {0.0035, 0.0} };



    // initial iteration 
    roots s_last = step(N, z_last);
    
    cout<<" ( "<<z.epsilon<<" ), i( "<<z.delta<<" ) ";
            
    for (int iter=0; iter<20; ++iter) {
        
        //roots z_next;
        roots s = step(N, z);
        
        s_last = s;
        z_last = z;
        z = s;
        
        cout<< "iter "<< iter <<" ";
        cout<<" ( "<<z.epsilon<<" ), i( "<<z.delta<<" ) ";
              
    }
/*
    cout<<endl<<" roots "<<get_roots(z)<<endl;
    
    cout<<endl<<" qn "<<dirtyQuantumJforRoots(N, get_roots(z))<<endl;
 */       
    return 0;    

}
