
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


struct roots {
    bool has_origin_root;
    bool has_singular_pair;
    vector<double> x;
    vector<double> y;    
    vector<double> kite_x;
    vector<double> kite_y;
    vector<double> bar_x;
    vector<double> bar_y;
};

vector< complex<double> > get_roots(const roots& in)
{
    
    vector< complex<double> > out;
    if (in.has_origin_root)  out.push_back(0.0);
    if (in.has_singular_pair) {
        out.push_back(0.5*I);
        out.push_back(-0.5*I);
    };
    for (int i = 0; i < in.x.size(); ++i) {
        out.push_back(in.x[i]);
        out.push_back(-in.x[i]);
    }
    for (int i = 0; i < in.y.size(); ++i) {
        out.push_back(I*in.y[i]);
        out.push_back(-I*in.y[i]);
    }
    for (int i = 0; i < in.kite_x.size(); ++i) {
        out.push_back(in.kite_x[i]);
        out.push_back(-in.kite_x[i]);
    }
    for (int i = 0; i < in.kite_y.size(); ++i) {
        out.push_back(I*in.kite_y[i]);
        out.push_back(-I*in.kite_y[i]);
    }
    for (int i = 0; i < in.bar_x.size(); ++i) {
        out.push_back(in.bar_x[i]+I*in.bar_y[i]);
        out.push_back(-in.bar_x[i]+I*in.bar_y[i]);
        out.push_back(in.bar_x[i]-I*in.bar_y[i]);
        out.push_back(-in.bar_x[i]-I*in.bar_y[i]);
    }
    
    return out;
}



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



roots step(const int chain_length, const vector<int>& bethe_2I_x, const vector<int>& bethe_2I_bar, const roots& in) {
    
    roots out = in;  
    // copy not necessary, would this be faster?
    // roots out = {vector<double> (in.x.size()), vector<double> (in.y.size())};
    
    // imaginary pairs   
    for (int i=0; i< in.y.size(); ++i) {
    
        double scatter = 0.0;
        for (int j=0; j< in.y.size(); ++j) {
            // this sum excludes the self-scattering term of the string, which is taken to the LHS  
            if (j==i) continue;
            scatter += atih(in.y[i] - in.y[j]) + atih(in.y[i] + in.y[j]);
        }
        for (int j=0; j< in.x.size(); ++j) {
             // scattering with real pairs
            scatter += zeta(in.x[j], 1.0+in.y[i]) - zeta(in.x[j], 1.0-in.y[i]);
            //atan (iy-x) + atan(iy+x) = atan(x+iy) - atan(x-iy) = i*(ln (sqrt( x^2 + (1+y)^2 )) - ln (sqrt(x^2 + (1-y)^2)) )
        }
        for (int j=0; j< in.kite_x.size(); ++j) {
             // scattering with kites
            scatter += zeta(in.kite_x[j], 1.0+in.y[i]) - zeta(in.kite_x[j], 1.0-in.y[i]);
            scatter += atih(in.y[i] - in.kite_y[j]) + atih(in.y[i] + in.kite_y[j]);
        }
        for (int j=0; j< in.bar_x.size(); ++j) {
            // scattering with bars
            scatter += zeta(in.bar_x[j], 1.0+in.y[i]+in.bar_y[j]) - zeta(in.bar_x[j], 1.0-in.y[i]-in.bar_y[j]);
            scatter += zeta(in.bar_x[j], 1.0+in.y[i]-in.bar_y[j]) - zeta(in.bar_x[j], 1.0-in.y[i]+in.bar_y[j]);
        }
        // scatter off origin root
        if (in.has_origin_root) scatter += atih(in.y[i]);
        // scatter off singular pair
        if (in.has_singular_pair) scatter += atih(in.y[i]-0.5)+atih(in.y[i]+0.5);
        
        // LHS == (N-1) atih (2y) == kinetic term + self scattering term.
        // this assumes wide pairs
        out.y[i] = 0.5/tanh( scatter / double(chain_length-1)  ) ; 
    }
    
    // real pairs
    for (int i=0; i< in.x.size(); ++i) {
        double scatter = 0.0;
        for (int j=0; j< in.y.size(); ++j) {;
            // scattering with imaginary pairs
            scatter += xi(in.x[i], 1.0+in.y[j]) + xi(in.x[i], 1.0-in.y[j]);
        }
        for (int j=0; j< in.x.size(); ++j) {
            // scattering with real pairs
            // this excludes self-scattering of the real pair, taken to the LHS 
            if (i==j) continue;
            scatter += atan(in.x[i]-in.x[j]) + atan(in.x[i]+in.x[j]);
        }
        for (int j=0; j< in.kite_x.size(); ++j) {
            // scattering with kites
            scatter += xi(in.x[i], 1.0+in.kite_y[j]) + xi(in.x[i], 1.0-in.kite_y[j]);
            scatter += atan(in.x[i] - in.kite_x[j]) + atan(in.x[i] + in.kite_x[j]);
        }
        for (int j=0; j< in.bar_x.size(); ++j) {
            // scattering with bars
            scatter += xi(in.x[i]-in.bar_x[j], 1.0+in.bar_y[j]) + xi(in.x[i]-in.bar_x[j], 1.0-in.bar_y[j]);
            scatter += xi(in.x[i]+in.bar_x[j], 1.0+in.bar_y[j]) + xi(in.x[i]+in.bar_x[j], 1.0-in.bar_y[j]);
        }
        // scatter off origin root
        if (in.has_origin_root) scatter += atan(in.x[i]);
        // scatter off singular pair
        if (in.has_singular_pair) scatter += xi(in.x[i], 1.5) + xi(in.x[i], 0.5);
                               
        // LHS == (N-1) atan (2x) == kinetic term + self-scattering 
        out.x[i] = 0.5*tan( (0.5*PI*bethe_2I_x[i] + scatter) / double(chain_length-1)  ) ; 
    }
    
    // kites
    for (int i=0; i< in.kite_x.size(); ++i) {
        complex<double> scatter = 0.0;
        for (int j=0; j< in.y.size(); ++j) {;
            // scattering with imaginary pairs
            scatter += xi(in.kite_x[i], 1.0+in.y[j]) + xi(in.kite_x[i], 1.0-in.y[j]);
            scatter += I*( atih(in.kite_y[i] - in.y[j]) + atih(in.kite_y[i] + in.y[j]) );
        }
        for (int j=0; j< in.x.size(); ++j) {
            // scattering with real pairs
            scatter += atan(in.kite_x[i]-in.x[j]) + atan(in.kite_x[i]+in.x[j]);
            scatter += I*( zeta(in.x[j], 1.0+in.kite_y[i]) - zeta(in.x[j], 1.0-in.kite_y[i]) );
        }
        for (int j=0; j< in.kite_x.size(); ++j) {
            // scattering with kites, excluding self-scattering
            if (i==j) continue;
            scatter += xi(in.kite_x[i], 1.0+in.kite_y[j]) + xi(in.kite_x[i], 1.0-in.kite_y[j]);
            scatter += atan(in.kite_x[i] - in.kite_x[j]) + atan(in.kite_x[i] + in.kite_x[j]);
            
            scatter += I*( zeta(in.kite_x[j], 1.0+in.kite_y[i]) - zeta(in.kite_x[j], 1.0-in.kite_y[i]) );
            scatter += I*( atih(in.kite_y[i] - in.kite_y[j]) + atih(in.kite_y[i] + in.kite_y[j]) );
        }
        for (int j=0; j< in.bar_x.size(); ++j) {
            // scattering with bars
            scatter += xi(in.kite_x[i]-in.bar_x[j], 1.0+in.bar_y[j]) + xi(in.kite_x[i]-in.bar_x[j], 1.0-in.bar_y[j]);
            scatter += xi(in.kite_x[i]+in.bar_x[j], 1.0+in.bar_y[j]) + xi(in.kite_x[i]+in.bar_x[j], 1.0-in.bar_y[j]);
            
            scatter += I*( zeta(in.bar_x[j], 1.0+in.kite_y[i]+in.bar_y[j]) - zeta(in.bar_x[j], 1.0-in.kite_y[i]-in.bar_y[j]));
            scatter += I*( zeta(in.bar_x[j], 1.0+in.kite_y[i]-in.bar_y[j]) - zeta(in.bar_x[j], 1.0-in.kite_y[i]+in.bar_y[j]));
        }
        
        // scatter off origin root
        if (in.has_origin_root) 
            scatter += atan(in.kite_x[i]) + I*atih(in.kite_y[i]);
        // scatter off singular pair
        if (in.has_singular_pair) { 
            scatter += xi(in.kite_x[i], 1.5) + xi(in.kite_x[i], 0.5);
            scatter += I*( atih(in.kite_y[i] - 0.5) + atih(in.kite_y[i] + 0.5) );
        }
        
        // sum of 2* positive-root quantum numbers - immaterial because of tan, as long as it's an even multiple of pi
        //scatter += PI*(double(chain_length)-2);
                            
        complex<double> kite_z = tan( 0.5*double(chain_length-1) * (atan(2.0*in.kite_x[i]) + atan(2.0*I*(in.kite_y[i]))) - 0.5*scatter);
                               
        out.kite_x[i] = real(kite_z);
        out.kite_y[i] = imag(kite_z);
         
    }
    
    // bars
    for (int i=0; i< in.bar_x.size(); ++i) {
        complex<double> scatter = 0.0;
        complex<double> z_i = in.bar_x[i] + I*in.bar_y[i];
        for (int j=0; j< in.y.size(); ++j) {;
            // scattering with imaginary pairs
            scatter += atan(z_i - I*in.y[j]) + atan(z_i + I*in.y[j]);
        }
        for (int j=0; j< in.x.size(); ++j) {
            // scattering with real pairs
            // this excludes self-scattering of the real pair, taken to the LHS 
            scatter += atan(z_i - in.x[j]) + atan(z_i + in.x[j]);
        }
        for (int j=0; j< in.kite_x.size(); ++j) {
            // scattering with kites
            scatter += atan(z_i - I*in.kite_y[j]) + atan(z_i + I*in.kite_y[j]);
            scatter += atan(z_i - I*in.kite_x[j]) + atan(z_i + I*in.kite_x[j]);
        }
        for (int j=0; j< in.bar_x.size(); ++j) {
            // scattering with bars
            if (i==j) continue;
            scatter += atan(z_i - (in.bar_x[j] + I*in.bar_y[j])) + atan(z_i - (in.bar_x[j] - I*in.bar_y[j]));
            scatter += atan(z_i + (in.bar_x[j] + I*in.bar_y[j])) + atan(z_i + (in.bar_x[j] - I*in.bar_y[j]));
        }
        // scatter off origin root
        if (in.has_origin_root) {
            scatter += atan(z_i);
        }
        // scatter off singular pair
        if (in.has_singular_pair) {
            scatter += atan(z_i - I*0.5) + atan(z_i + I*0.5);
        }
                               
        // LHS == (N-1) atan (2x) == kinetic term + self-scattering 
        complex<double> bar_z = 0.5*tan( (0.5*PI*bethe_2I_bar[i] + atan(2.0*I*in.bar_y[i]) + atan(2.0*in.bar_x[i]) + scatter) 
                                    / double(chain_length-1)  ) ; 
        out.bar_x[i] = real(bar_z);
        out.bar_y[i] = imag(bar_z);

    }
    
    return out;
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

    
    // NOTE: ys should be strictly more than 1 apart
    
/*
    //OK
    //iter 23  (  ), i( 1.00001480038088 2.16376688302984 )  (  ), i(  )  conv 2.70076219416421e-24
    //#6#0: 2I [;;;; +0 ;]  2J [-8 -14 -2 -10 -6] sums [+0;] sum +0 (-3.1791694e-15,+2.1637669) (-1.8111943e-15,+1.0000148) (-1.8106887e-15,+0) (-1.8111943e-15,-1.0000148) (-3.1791694e-15,-2.1637669)
    int N = 10;
    vector<int> bethe_2I_x = {};
    roots z_last = { true, false, {}, {1.1, 2.2} };
    roots z = { true, false, {}, {1.2, 2.3} };
*/

/*
    //OK
    //iter 10  ( 0.216201306054027 0.598086994989464 ), i(  )  (  ), i(  )  conv 1.34717932379393e-23
    //#0#0: 2I [ +0 +2 -2 +4 -4 ;]  2J [+0 +2 -2 +4 -4] sums [+0;+2;-2;+4;-4;] sum +0 (+0,+0) (+0.21620131,+0) (-0.21620131,+0) (+0.59808699,+0) (-0.59808699,+0)
    int N = 10;
    vector<int> bethe_2I_x = {2, 4};
    roots z_last = { true, false, {0.01, 1.1}, {} };
    roots z = { true, false, {0.009, 1.2}, {} };
*/
/*
    // single kite
    // OK
    //#3#0: 2I [ +0 ;; +0 ;]  2J [-1 +9 +1 -9] sums [-1;+1;] sum +0 (+0.01414495,+0) (+0,+1.0035739) (-0.01414495,+0) (+0,-1.0035739)
    //iter 8  (  ), i(  )  ( 0.0141449501592471 ), i( 1.00357387588898 )  conv 2.27585135999056e-21
    int N = 10;
    vector<int> bethe_2I_x = {};
    vector<int> bethe_2I_bar = {};

    roots z_last = { false, false, {}, {}, {0.014},  {1.0034} };
    roots z = { false, false, {}, {}, {0.0141}, {1.0035} };
 */   
 /*
    // single kite, alternative
    // this does not work.
    //#3#0: 2I [ +0 ;; +0 ;]  2J [-1 +9 +1 -9] sums [-1;+1;] sum +0 (+0.01414495,+0) (+0,+1.0035739) (-0.01414495,+0) (+0,-1.0035739)
    //iter 8  (  ), i(  )  ( 0.0141449501592471 ), i( 1.00357387588898 )  conv 2.27585135999056e-21
    int N = 10;
    vector<int> bethe_2I_x = {-1};
    vector<int> bethe_2I_bar = {};

    roots z_last = { false, false,  {0.014},  {1.0034} };
    roots z = { false, false,  {0.0141}, {1.0035} };
   */ 
    
    

/*
    //OK
    //iter 5  ( 0.23612398353546 ), i(  )  (  ), i(  )  conv 1.76019828615352e-28
    //#1#0: 2I [ +0 +2 -2 ; +0 ;]  2J [+0 +2 -2 +4 +6] sums [+0;+2;-2;+10;] sum +10 (+4.2284582e-21,+0) (+0.23612398,+0) (-0.23612398,+0) (-2.1684043e-20,+0.5) (-2.1684043e-20,-0.5)
    int N = 10;
    vector<int> bethe_2I_x = {2};
    roots z_last = { true, true, {1.1}, {} };
    roots z = { true, true, {1.2}, {} };
*/

    // OK
    //#3#0: 2I [ +1 -1 ;; +0 ;]  2J [+0 +0 +8 +2 +10] sums [+0;+0;+0;] sum +0 (+0.11917111,+0) (-0.11917111,+0) (+1.2598685e-16,+1.0446068) (+5.534098e-17,+0) (+1.2598685e-16,-1.0446068)
    //iter 100  ( 0.119171104825923 ), i( 1.04460675205512 )  (  ), i(  )  conv 2.56131024799912e-22
    int N = 10;
    vector<int> bethe_2I_x = {0};
    vector<int> bethe_2I_bar = {};

    roots z_last = { true, false, {1.1}, {1.1} };
    roots z = { true, false, {1.2}, {1.2} };

/*
    // kite + singular pair
    // roots (0,0.5) (-0,-0.5) (-0.0182952672265357,0) (0.0182952672265357,0) (0,1.00656301253971) (-0,-1.00656301253971)
    // qn (2.5,0) (3.5,0) (0.500000000677972,0) (-0.500000000677972,0) (5,3.0734837818768e-10) (-5,-3.0734837818768e-10)
    // OK
    int N = 12;
    vector<int> bethe_2I_x = {};
    vector<int> bethe_2I_bar = {};
    roots z_last = { false, true, {}, {}, {1.014},  {1.0034} };
    roots z = { false, true, {}, {}, {1.0141}, {1.0035} };
*/
/*
    // 8-string
    // OK 
    // roots (0,0.5) (-0,-0.5) (0,1.50000218850042) (-0,-1.50000218850042) (0,2.52077526622762) (-0,-2.52077526622762) (0,4.49648402696482) (-0,-4.49648402696482)
    // qn (3.5,0) (4.5,0) (6.5,7.47447136404745e-12) (-6.5,-7.47465689640387e-12) (5.5,2.14258772906719e-11) (-5.5,-2.14259479696648e-11) (4.5,6.9150559699486e-14) (-4.5,-6.91682294477163e-14)
    // second qn should be -3.5
    
    int N = 16;
    vector<int> bethe_2I_x = {};
    vector<int> bethe_2I_bar = {};
    roots z_last = { false, true, {}, {1.6, 2.7, 3.8}  };
    roots z = { false, true, {}, {1.61, 2.71, 3.81} };
*/

/*
    int N = 16;
    vector<int> bethe_2I_x = {};
    vector<int> bethe_2I_bar = {9};
    roots z_last = { false, false, {}, {}, {0.5}, {2.15}, {0.1}, {1.25}  };
    roots z = { false, false, {}, {}, {0.6}, {2.2}, {0.11}, {1.3} };
*/
/*
    int N = 12;
    vector<int> bethe_2I_x = {};
    vector<int> bethe_2I_bar = {};
    roots z_last = { false, false, {}, { 1.1, 2.2, 3.3}  };
    roots z = { false, false, {}, { 1.11, 2.21, 3.31} };
*/

    // initial iteration 
    roots s_last = step(N, bethe_2I_x, bethe_2I_bar, z_last);
    
    // enforce wide strings
    // minimum: i, 2i, etc for odd numbers of pairs on the imaginary axis
    // 0.5i, 1.5i, etc for even numbers
     
    double shift = (z.has_singular_pair) ? 1.5 : 1.0;
            
    for (int iter=0; iter<200; ++iter) {
        double convergence = 0;
        
        roots z_next;
        z_next.has_origin_root = z.has_origin_root;
        z_next.has_singular_pair = z.has_singular_pair;

        roots s = step(N, bethe_2I_x, bethe_2I_bar, z);
        
        z_next.y = secant_step(z_last.y, z.y, s_last.y, s.y, convergence);
        z_next.x = secant_step(z_last.x, z.x, s_last.x, s.x, convergence);
        
        z_next.kite_x = secant_step(z_last.kite_x, z.kite_x, s_last.kite_x, s.kite_x, convergence);
        z_next.kite_y = secant_step(z_last.kite_y, z.kite_y, s_last.kite_y, s.kite_y, convergence);
        
        z_next.bar_x = secant_step(z_last.bar_x, z.bar_x, s_last.bar_x, s.bar_x, convergence);
        z_next.bar_y = secant_step(z_last.bar_y, z.bar_y, s_last.bar_y, s.bar_y, convergence);
        
        
        for (int i = 0; i < z_next.y.size(); ++i) {
            // enforce wide strings
            // minimum: i, 2i, etc for odd numbers of pairs on the imaginary axis
            // 0.5i, 1.5i, etc for even numbers 
            if (z_next.y[i] <= shift+i) {
                z_next.y[i]= 0.5*(shift + i + z.y[i]);
                // indicate no convergence
                convergence += 1; 
            }
        }


        for (int i = 0; i < z_next.bar_y.size(); ++i) {
            if (z_next.bar_y[i] < 1.0) {
                z_next.bar_y[i]= 0.5*(1.0 + z.bar_y[i]);
                convergence += 1; 
            }
        }
        for (int i = 0; i < z_next.bar_x.size(); ++i) {
            if (z_next.bar_x[i] <= 0.0) {
                z_next.bar_x[i]= 0.5*(z.bar_x[i]);
                convergence += 1; 
            }
        }
    
        
        
        s_last = s;
        z_last = z;
        z = z_next;
        
        cout<< "iter "<< iter <<" ";
        cout<<" ( "<<z_next.x<<" ), i( "<<z_next.y<<" ) ";
        cout<<" ( "<<z_next.kite_x<<" ), i( "<<z_next.kite_y<<" ) ";
        cout<<" ( "<<z_next.bar_x<<" )+i( "<<z_next.bar_y<<" ) "; 
        cout<<" conv "<<convergence<<endl;
        if (convergence < 1e-20) break;        
    }

    cout<<endl<<" roots "<<get_roots(z)<<endl;
    
    cout<<endl<<" qn "<<dirtyQuantumJforRoots(N, get_roots(z))<<endl;
        
    return 0;    

}
