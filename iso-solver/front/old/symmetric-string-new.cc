#include "symmetric-string.h"

using namespace std;


/** class String **/


String::String(const int big_n, const int string_length, const int sum_jx2, const double lambda)
      : big_n_(big_n), 
        string_length_(string_length), sum_jx2_(sum_jx2),  number_roots_odd_(false),
        lambda_(lambda), epsilon_(string_length), delta_(string_length), 
        new_lambda_(0.), new_epsilon_(string_length), new_delta_(string_length), 
        hold_(string_length, false),
        iterations_(0) 
{
    setPolicy();
}

void String::setPolicy(const double initial_deviation, const int steps_no_deviation, const double damping_delta)
{
    initial_deviation_ = initial_deviation;
    steps_no_deviation_ = steps_no_deviation;
    damping_delta_ = damping_delta;
}
    

vector< complex<double> > String::getComplexQuartets() const
{
    vector< complex<double> > result(string_length_/2);
    // run through half string, excluding odd strings' centres which are real.
    for (int a=0; a<string_length_/2; ++a) {
        result[a] = lambda_ + epsilon_[a] + I*(delta_[a] + 0.5*(string_length_-1) - a);
    }
    return result;
}

vector< double > String::getRealPairs() const
{
    // no real root for even string
    if (string_length_%2) return { lambda_ };
    else return { };
}

    
vector< complex<double> > String::getRoots() const
{
    vector< complex<double> > roots (string_length_*2);
    for (int a=0; a<string_length_; ++a) {
        // run through string, this covers conjugates.
        roots[a] = lambda_ + epsilon_[a] + I*(delta_[a] + 0.5*(string_length_-1) - a);
        // opposite string
        roots[string_length_*2-a-1] = -roots[a];
    }
    return roots;
}

int String::size() const
{
    return string_length_*2;
}


long double String::thetaOthers(const double x, const double y, 
                            bool has_origin, const vector<double>& other_x, 
                            const vector<double>& other_y, const vector<complex<double> >& other_z) 
                            
{
	long double theta = 0.;
	
    if (has_origin) {
        theta +=  xi(x, 1.+y) + xi(x, 1.-y);
    }
    
    for (int k=0; k<other_x.size(); ++k) {
        double x_k = other_x[k];
	    theta +=  xi( x-x_k, 1.+y )  +  xi( x-x_k, 1.-y );
	    theta +=  xi( x+x_k, 1.+y )  +  xi( x+x_k, 1.-y );
    }
    
    for (int k=0; k<other_y.size(); ++k) {
        double y_k = other_y[k];
	    theta +=  xi( x, 1.+(y-y_k) )  +  xi( x, 1.-(y-y_k) );
	    theta +=  xi( x, 1.+(y+y_k) )  +  xi( x, 1.-(y+y_k) );
    }
    
    for (int k=0; k<other_z.size(); ++k) {
        double x_k = real(other_z[k]);
        double y_k = imag(other_z[k]);
	    theta +=  xi( x-x_k, 1.+(y-y_k) )  +  xi( x-x_k, 1.-(y-y_k) );
	    theta +=  xi( x+x_k, 1.+(y+y_k) )  +  xi( x+x_k, 1.-(y+y_k) );
	    theta +=  xi( x-x_k, 1.+(y+y_k) )  +  xi( x-x_k, 1.-(y+y_k) );
	    theta +=  xi( x+x_k, 1.+(y-y_k) )  +  xi( x+x_k, 1.-(y-y_k) );
    }

    return theta;
}


long double String::logRSqOthers(const double x, const double y, 
                            bool has_origin, const vector<double>& other_x, 
                            const vector<double>& other_y, const vector<complex<double> >& other_z) 
{
	long double log_rsq = 0.;
	
    if (has_origin) {
        log_rsq +=  log( sq(x) + sq(1.+y) ) - log( sq(x) + sq(1.-y) );
    }
    
    for (int k=0; k<other_x.size(); ++k) {
        double x_k = other_x[k];
		log_rsq +=  log( sq(x-x_k) + sq(1.+y) ) - log( sq(x-x_k) + sq(1.-y) );
		log_rsq +=  log( sq(x+x_k) + sq(1.+y) ) - log( sq(x+x_k) + sq(1.-y) );
    }
    
    for (int k=0; k<other_y.size(); ++k) {
        double y_k = other_y[k];
		log_rsq +=  log( sq(x) + sq(1.0 + (y-y_k)) ) - log( sq(x) + sq(1.0 - (y-y_k)) );
		log_rsq +=  log( sq(x) + sq(1.0 + (y+y_k)) ) - log( sq(x) + sq(1.0 - (y+y_k)) );
    }
    
    for (int k=0; k<other_z.size(); ++k) {
        double x_k = real(other_z[k]);
        double y_k = imag(other_z[k]);
		log_rsq +=  log( sq(x-x_k) + sq(1.0 + (y-y_k)) ) - log( sq(x-x_k) + sq(1.0 - (y-y_k)) );
		log_rsq +=  log( sq(x+x_k) + sq(1.0 + (y+y_k)) ) - log( sq(x+x_k) + sq(1.0 - (y+y_k)) );
 		log_rsq +=  log( sq(x-x_k) + sq(1.0 + (y+y_k)) ) - log( sq(x-x_k) + sq(1.0 - (y+y_k)) );
 		log_rsq +=  log( sq(x+x_k) + sq(1.0 + (y-y_k)) ) - log( sq(x+x_k) + sq(1.0 - (y-y_k)) );
    }

    return log_rsq;
}




/** calculate a new value for the deviance (delta) and aberration (epsilon) **/
bool String::stepDeviation (
                const vector<SymRoots*>& all, const int alpha, 
                vector<double>& new_delta, vector<double>& new_epsilon)
{
	// real rapidity, no deviations.
	if (1==string_length_) 		{
	    new_delta[0] = new_epsilon[0] = { 0. };
		return true;
	}

	long double theta = 0.0;
	long double log_rsq = 0.0;

	const int a_mid = (string_length_-1)/2;
	long double del_term[string_length_/2];
	long double eps_term[string_length_/2];
	
	for (int a=0; a<string_length_/2; ++a) {
		 
		const double x = lambda_ + epsilon_[a];
	    const double y = 0.5*(string_length_-1) - a + delta_[a];
	    
		// string self scattering minus dangerous terms		
		for (int b=0; b < string_length_; ++b) {
			double x_sum = 2.*lambda_ + epsilon_[a] + epsilon_[b];
			double y_sum = (string_length_-1.) - (a+b) + delta_[a] + delta_[b];
	        double x_diff = epsilon_[a] - epsilon_[b];
	        double y_diff = (delta_[a]-delta_[b]) - (a-b);
			if (b!=a) {    
			    if (a-b!=1)  theta -=  xi( x_diff, 1.+y_diff );
			    if (a-b!=-1) theta -=  xi( x_diff, 1.-y_diff );
			    if (a-b!=1)  log_rsq +=  log( sq(x_diff) + sq(1.+y_diff) );
			    if (a-b!=-1) log_rsq += -log( sq(x_diff) + sq(1.-y_diff) );
		        // opposite string, excluding opposite root absorbed in kinetic term
			    theta -=  xi( x_sum, 1.-y_sum ) + xi( x_sum, 1.+y_sum );
			    log_rsq +=  log( sq(x_sum) + sq(1.+y_sum) ) - log( sq(x_sum) + sq(1.-y_sum) );
		    }
		}		
	    
	    for (int beta=0; beta<all.size(); ++beta) {
	        // exclude self scattering
	        if (beta==alpha) continue;
	    	// scattering with others
	    	const bool other_o = all[beta]->hasOrigin();
	    	const vector<double> other_x = all[beta]->getRealPairs();
	    	const vector<double> other_y = all[beta]->getImagPairs();
	    	const vector<complex<double> > other_z = all[beta]->getComplexQuartets();
            theta -= thetaOthers( x, y, other_o, other_x, other_y, other_z );
            log_rsq += logRSqOthers( x, y, other_o, other_x, other_y, other_z );
        }
            
		// contribution M mod 2 from sum of conjugate quantum numbers
		theta -= PI*number_roots_odd_; 
		
		// kinetic theta, absorbing opposite root
		theta += (big_n_-1.) * ( xi(2.*x, 1.+2.*y) + xi(2.*x, 1.-2.*y) );
		log_rsq += (big_n_-1.) * ( log(sq(x) + sq(y-0.5)) - log(sq(x) + sq(y+0.5)) );

		// theta is now sum b=1..a theta[b]
		// log_r is now sum b=1..a log_r[b]
		del_term[a] = -cos(theta)*exp(0.5*log_rsq);
		eps_term[a] = sin(theta)*exp(0.5*log_rsq);
	}	
		
    // outputs
    new_delta.assign(string_length_, 0.);
    new_epsilon.assign(string_length_, 0.);

    // add differences
    for (int a=a_mid-1; a>=0; --a) {   
	    new_delta[a] += new_delta[a+1] + del_term[a];
        new_epsilon[a] += new_epsilon[a+1] + eps_term[a];
    }

	//new_epsilon[a_mid] = 0.0;
	//new_delta[a_mid] = 0.0;
	if (0==string_length_%2) {
		// even n_j: apply the inner-sign term.
        bool inner_narrow = innerPairIsNarrow(sum_jx2_, string_length_, number_roots_odd_);
		for (int a=0; a<a_mid+1; ++a) {
			// log_r should be sum b=1...n/2 log_r[b] now
			new_delta[a] += 0.5*exp(0.5*log_rsq) * (inner_narrow?-1.:1.); 
		}
	}

	// set the conjugates
	for (int a=a_mid+1; a<string_length_; ++a) {
		new_delta[a] = - new_delta[string_length_-a-1];
		new_epsilon[a] = new_epsilon[string_length_-a-1];
	}

	return true;
}




/** calculate a new value for the string centre (lambda) **/
double String::stepRapidity (const vector<SymRoots*>& all, const int alpha)
{
	long double scatter = 0.0;		
	for (int a=0; a < string_length_; ++a) {

        const double x = lambda_ + epsilon_[a];
        const double y = 0.5*(string_length_-1) - a + delta_[a];
            
        // scattering with opposite strings
	    for (int b=0; b < string_length_; ++b) {
            // direct opposite - absorb into kinetic term
            if (b+a==string_length_-1) continue;
 
            const double x_b = - lambda_ - epsilon_[b];
	        const double y_b = 0.5*(string_length_-1) - b + delta_[b];
 			scatter += xi( x-x_b, 1.-(y-y_b) )  +  xi( x-x_b, 1.+(y-y_b) );
		}

        for (int beta=0; beta<all.size(); ++beta) {
	        // exclude self scattering
	        if (beta==alpha) continue;
	    	// scattering with others
	    	const bool other_o = all[beta]->hasOrigin();
	    	const vector<double> other_x = all[beta]->getRealPairs();
	    	const vector<double> other_y = all[beta]->getImagPairs();
	    	const vector<complex<double> > other_z = all[beta]->getComplexQuartets();
            scatter += thetaOthers( x, y, other_o, other_x, other_y, other_z );
        }
        
	}
	
	// 2\pi \sum_{a=1}^n J_a + \theta\sub{other}^a = \sum_{a=1}^n \theta\sub{kin}^a 
    // absorbed opposite root in N-1 for faster convergence
    const long double phase = 0.5*(PI*sum_jx2_ + scatter) / (big_n_-1.);
    const int a_mid = (string_length_-1)/2;

	long double sum_kin = 0.0;
	// term for even strings
	if (0==string_length_%2)   sum_kin += xi( lambda_, -delta_[a_mid]); 
	
    for (int a=0; a< a_mid; ++a) {
		sum_kin += 	  xi( lambda_ + epsilon_[a+1], 0.5*string_length_ - (a+1) + delta_[a+1]  )
					+ xi( lambda_ + epsilon_[a],  -0.5*string_length_ + (a+1) - delta_[a]  );
	}

    return -epsilon_[0] + ( 0.5*string_length_ + delta_[0] )*tan( phase-sum_kin ) ;
}



bool String::iterate(const vector<SymRoots*>& all, const int alpha, const double last_convergence, double& convergence)
{
	new_lambda_ = stepRapidity(all, alpha);
	if (!finite(new_lambda_)) return false;
	 
    if (iterations_ > steps_no_deviation_) {
	    if (!stepDeviation(all, alpha, new_delta_, new_epsilon_)) return false;
        for (int a=0; a<string_length_; ++a) 
	        new_delta_[a] = dampedStep(delta_[a], new_delta_[a], damping_delta_);
    }
    
    // enforce holds
    for (int a=0; a<string_length_;++a) 
        if (hold_[a]) {
		    new_delta_[a] = delta_[a];
		    new_epsilon_[a] = epsilon_[a];
	    }
    
    // calculate convergence
	convergence += sq(lambda_-new_lambda_);
	for (int a=0; a < string_length_; ++a)
		convergence += sq(new_delta_[a] - delta_[a]) + sq(new_epsilon_[a] - epsilon_[a]);

    ++iterations_;
	// report success	
	return true;
}



void String::refresh()
{
	lambda_ = new_lambda_;
	delta_ = new_delta_;
	epsilon_ = new_epsilon_;
    return;
}

	
bool String::initiate(const vector<SymRoots*>& all, const int alpha)
{
    number_roots_odd_ = false;
    for (int beta=0; beta<all.size();++beta) { 
        if (all[beta]->hasOrigin()) {
            number_roots_odd_=true;
            break;
        } 
    }
    
    if (0.==lambda_) {
        if (0==sum_jx2_)
            // initial value must not be zero 
            // I hope this will stay out of the way if the int*Pi or half-int*Pi for the others.
            lambda_ = 1e-2/big_n_;
        else
            // non-interacting solution
            // this would need bethe-takahashi qn, not sum_jx2_!
            // in current form, leads to negative starting points which lead to non-convergence. 
            //lambda_ = 0.5*string_length_*tan(0.5*sum_jx2_*PI/N_);    
            
            //lambda_ = 0.25*sum_jx2_*PI/(N_-1.);    
            // not a good guess, but at least different for string with different qns. 
            lambda_ = abs(0.5*sum_jx2_/double(string_length_*big_n_));    
    }
    
    const int a_mid = (string_length_-1)/2;
    for(int a=0; a<a_mid; ++a) {
        epsilon_[a] = epsilon_[string_length_-1-a] = initial_deviation_;
        delta_[a] = initial_deviation_;
        delta_[string_length_-1-a] = -delta_[a];
    }
    epsilon_[a_mid] = 0.;
    	
    // for odd strings, set central deviations.
    // for even strings, ensure the deviation of the inner pair has the correct sign
    if (string_length_%2) {
        delta_[a_mid] = 0.;
    }
    else {
        if (innerPairIsNarrow(sum_jx2_, string_length_, number_roots_odd_)) {
            delta_[a_mid] = -initial_deviation_;
            delta_[a_mid+1] = -delta_[a_mid];
        } 
    }
    
    return true;
}



