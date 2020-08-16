
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




/* calculate the quantum numbers (J) as they are in the complex (non-Takahashi) Bethe equations, from the solution */
vector< complex<double> > cleanerJ(const double lambda, const double epsp, const double epsm, const double del,
                const double big_n)
{
    
 	vector< complex<double> > bethe_i (4);
 	
	//for (int alpha=0; alpha < for_roots.size(); ++alpha) {
	bethe_i[0] = atan(2.0*(lambda+I+I*del) ) * big_n - atan( I+I*del - epsp ) - atan( I+I*del + epsm ) - atan (2.*I+2.*I*del);
cerr<<atan(2.0*(lambda+I+I*del) ) * big_n/PI<<endl;
	cerr<<- atan( I+I*del - epsp )/PI<<endl;
	cerr<< - atan( I+I*del + epsm )/PI<<endl;
	cerr<< - atan (2.*I+2.*I*del)/PI<<endl;
	
	bethe_i[1] = atan(2.0*(lambda+epsp) ) * big_n - xi( epsp,  2.+del ) - xi (epsp, -del) - atan (epsp+epsm);
	
	bethe_i[2] = atan(2.0*(lambda-epsm) ) * big_n - xi( -epsm,  2.+del ) - xi (-epsm, -del) - atan (-epsp-epsm);
	
	bethe_i[3] = atan(2.0*(lambda-I-I*del) ) * big_n - atan( -I-I*del - epsp ) - atan( -I-I*del + epsm ) - atan (-2.*I-2.*I*del);
	   	
		
    bethe_i[0] /= PI;	
    bethe_i[1] /= PI;	
    bethe_i[2] /= PI;	
    bethe_i[3] /= PI;	
	
	return bethe_i;
}

/** calculate a new value for the string centre (lambda) **/
double stepRapidity (const double lambda_, const double epsp, const double epsm, const double del,
                const double big_n, double& new_lambda)  
{
    
	long double scatter = 0.0;		
    scatter += xi( epsm, -del ); //if eps=1e-10, del=1e-5 -> pi; if eps=1e-5, del=1e-10 -> pi/2  
    scatter += xi( epsm, 2.+del ); // 0
    scatter += atan(epsp+epsm); // 0
    
    double sum_jx2_ = -3.; //-3.;	    
	// -5 works with epsp=epsm >> del
	
	
	// 2\pi \sum_{a=1}^n J_a + \theta\sub{other}^a = \sum_{a=1}^n \theta\sub{kin}^a 
    long double phase = (0.5*PI*sum_jx2_ + scatter) / big_n;
    //phase = (-2 pi + 0.5pi [if del>eps])/12 

	long double sum_kin = 	 // atan( 2.*(lambda_ + epsp))
					+ xi( 2.*lambda_,  3.+2.*del  ) 
					+ xi( 2.*lambda_,  -1.-2.*del  )
					;
    // atan(tan(pi/3)/3) = pi/6 
    // + atan(tan(pi/3)/-1) 

	new_lambda =  0.5*tan( phase-sum_kin ) - epsp;
	//new_lambda =  ( 0.5*3. + del)*tan( phase-sum_kin ) ;
		
		
//    lambda = -0.5*tan(PI/3.);
// is a solution for small eps, delta only if del<eps


//cerr<<sum_kin/PI<<endl;	
//cerr<<"nl "<<new_lambda<<endl;
	return true;
}





/** calculate a new value for the deviance (delta) and aberration (epsilon) **/
bool stepDeviation (const double lambda_, const double epsp, const double epsm, const double del,
                const double big_n,  
                double& new_del, double& new_epsp, double& new_epsm)
{
    const int string_length_=3;
    //    const double    lambda_ = -0.5*tan(PI/3.);
    


	long double theta = 0.0;
	long double log_r = 0.0;
    
	long double thetap = 0.0;
	long double log_rp = 0.0;
    
	long double thetam = 0.0;
	long double log_rm = 0.0;
    
	
    // 4-> 2 (0 1 2 3),   3 -> 1 (0 1 2)
	//for (int a=0; a<string_length_/2; ++a) 
	//int a=0;
	{
		
		  
		// string self scattering minus dangerous terms		
		
		// with epsp
		
        double x_diff = - epsp;
        double y_diff = 1. + del;
    
	    //if (a-b!=1)  
	        theta -=  xi( x_diff, 1.+y_diff );
	    //if (a-b!=-1) 
	    //    thetam -=  xi( x_diff, 1.-y_diff );

	    //if (a-b!=1)  
	        log_r +=  log( sq(x_diff) + sq(1.+y_diff) );
	    //if (a-b!=-1) 
	    //    log_rm += -log( sq(x_diff) + sq(1.-y_diff) );
    
    
        // with epsm
        
        x_diff = epsm;
        
	    //if (a-b!=1)  
	        theta -=  xi( x_diff, 1.+y_diff );
	    //if (a-b!=-1) 
	    //    thetap -=  xi( x_diff, 1.-y_diff );

	    //if (a-b!=1)  
	        log_r +=  log( sq(x_diff) + sq(1.+y_diff) );
	    //if (a-b!=-1) 
	    //    log_rp += -log( sq(x_diff) + sq(1.-y_diff) );
    
    
	    // with -del
	    x_diff = 0.;
	    y_diff = 2.;
	
	    //if (a-b!=1)  
	        theta -=  xi( x_diff, 1.+y_diff );
	    //if (a-b!=-1) 
	        theta -=  xi( x_diff, 1.-y_diff );

	    //if (a-b!=1)  
	        log_r +=  log( sq(x_diff) + sq(1.+y_diff) );
	    //if (a-b!=-1) 
	        log_r += -log( sq(x_diff) + sq(1.-y_diff) );
    
        
		// contribution M mod 2 from sum of conjugate quantum numbers
	
	    const int number_roots = 4;
    	theta -= PI*number_roots;  
		
		// kinetic theta
		//double x = lambda_;
	    //double y = 1. + del;
	  
		//theta += big_n * ( xi(2.*x, 1.+2.*y) + xi(2.*x, 1.-2.*y) );
		theta += big_n * ( xi(2.*lambda_, 3.+2.*del) + xi(2.*lambda_, -1.-2.*del) );
		
		//log_r += big_n * ( log(sq(x) + sq(y-0.5)) - log(sq(x) + sq(y+0.5)) );
        log_r -= big_n * ( log(sq(2.*lambda_) + sq(3.+2.*del)) - log(sq(2.*lambda_) + sq(-1.-2.*del)) );

		double c = cos(theta)*exp(0.5*log_r);
		double d = sin(theta)*exp(0.5*log_r);
		
        
        log_r = 0.;
            
        // quantum numbers
        //double Jd = -1.5;
        //double Je = -4.5;
        //double Jp = -1.5; //guess
    
        //double x = lambda_+epsp;
        //y = 0.;
        
        double thetas =0. ;
        theta = 0.;
        theta += big_n* atan(2.*(lambda_+epsp)) - 0.5*PI;
        thetas += big_n* atan(2.*(lambda_+epsp)) - 0.5*PI;
   
   
   //cerr<<theta/PI+4.<<" ";     
        theta -= atan(epsp+epsm);
        
   //cerr<<theta/PI+4.<<" ";     
        theta -= xi(epsp, del+2.);
        thetas -= xi(epsp, del+2.) + xi(epsp, -del);
            
   //cerr<<theta/PI+4.<<" ";     
        // now theta = xi(epsp, -del)    
        double a = tan(theta);
        
        //new_epsp = -lambda_+ 0.5*tan((1./big_n) * (-4.5*PI + xi(epsp, -del) +xi(epsp, del+2.) + atan(epsp+epsm)));
        
        //x = lambda_-epsm;
        //y = 0.;
        double big_theta = theta;
        
        theta = 0.;
        theta += big_n* atan(2.*(lambda_-epsm)) - 0.5*PI;
        thetas -= big_n* atan(2.*(lambda_-epsm)) - 0.5*PI;
            
   //cerr<<theta/PI+4.<<" ";     
        theta -= atan(-epsp-epsm);
   //cerr<<theta/PI+4.<<" ";     
        theta -= xi(-epsm, del+2.);
        thetas += xi(-epsm, del+2.) + xi(-epsp, -del);
        
   //cerr<<theta/PI+4.<<" ";
        // now theta = xi(-epsm, -del) 
        double b = tan(theta);
        
        //new_epsm = lambda_- 0.5*tan((1./big_n) * (-3.5*PI + xi(-epsm, -del) +xi(-epsm, del+2.) + atan(-epsp-epsm)));
        
        big_theta +=theta;
        cerr<<"tanTh "<<tan(big_theta)<<endl;
        cerr<<"a "<<a<<endl;
cerr<<"b "<<b<<endl;
cerr<<"ab "<<a*b<<endl;
cerr<<"c "<<c<<endl;
cerr<<"d "<<d<<endl;    

        //new_del = sqrt(abs(-c + new_epsp*new_epsm));

        //new_del = 0.5*(-a*a+sqrt(a*a*a*a-4.*(a*d-c)));   

        //new_del = sgn(b)*sqrt(c/(1.-a*b));
        //new_epsp = -a*new_del;
        //new_epsm = +b*new_del; // d/new_del + new_epsp;

        double e = tan(thetas);
cerr<<" epsp+epsm " <<e<<endl;
        new_del = sgn(d)*sqrt(((4.*c-e*e)-sqrt(sq(4.*c-e*e)-16.*d*d))/8.);
        
        
        cerr<<" epsm-epsp "<<d/new_del<<endl;
        new_epsm = 0.5*(e+d/new_del);
        new_epsp = 0.5*(e-d/new_del);

cerr<<"del "<<new_del<<endl;
cerr<<"eps+ "<<new_epsp<<endl;
cerr<<"eps- "<<new_epsm<<endl;

	}	

	
	// signal success	
	return true;
}



 
          
int main() {
    // 0.866025
double    lambda = -0.5*tan(PI/3.);

    double N=12.0;
   
    double epsp = 1e-20;//1e-5; 
    double epsm = -1e-6;//1e-5; 
    
    double del = 1e-8;
    
    
        


    for(int i=0; i<10; ++i) {
        
        /*
        double theta = PI*(1.+Je) - N*atan(2.*lambda+2.*eps) + atan(2.*eps) + xi(eps, 2.+del);
        double r = -N*0.5*(sqrt(4.*sq(lambda)+sq(3.+2.*del)) - sqrt(4.*sq(lambda)+sq(1.+2.*del))) + 2.*atih(2.+2.*del) + 0.5*sqrt(sq(eps)+sq(2.+del));
        double eps_n = r*sin(theta);
        double del_n = r*cos(theta);
  cerr<<r<<" "<<theta<<" "<<eps_n<<" "<<del_n<<endl;  
        */
        
        double new_del=0.;
        double new_epsp=0.;
        double new_epsm=0.;
        double new_lambda = lambda;
        //if (i<=4) 
        {
            stepRapidity (lambda, epsp, epsm, del, N,
                 new_lambda)
            ;
            
            lambda = new_lambda;
        }
        
        double convergence = sq(lambda-new_lambda);
        if (i>4) {
            stepDeviation (lambda, epsp, epsm, del, N,
                     new_del,  new_epsp,  new_epsm)
            ;    
            convergence += sq(epsp-new_epsp) + sq(epsm-new_epsm) + sq(del-new_del) ;

            epsp = new_epsp;
            epsm = new_epsm;
            
            del = new_del;
        
        }
        
        lambda = new_lambda;
        
        cout.precision(4);
        
        cout<<endl;
        cout<<" iter "<<i;
        cout<<" lcnv "<<convergence;
            
        vector<complex<double> > roots(4);
        roots[0] = lambda + I*(1.+del);
        roots[1] = lambda - epsm ;
        roots[2] = lambda + epsp ;
        roots[3] = conj(roots[0]);
        
        cout<<" berr "<<betheError(N, roots);
        cout<<endl;    
        
        cout.precision(16);
        cout<<"rts "<<roots<<endl;
        cout.precision(4);
        cout<<"qns  "<<dirtyQuantumJforRoots(N, roots)<<endl;
        cout<<cleanerJ(lambda, epsp, epsm, del, N)<<endl;
    }    

    
    return 0;
}




