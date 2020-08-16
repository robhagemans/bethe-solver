#include "symmetric-string.h"

using namespace std;


/** class String **/


String::String(const int big_n, const int string_length, const int sum_jx2)// /*const int number_roots_total, const bool inner_narrow*/
      : big_n_(big_n), 
        string_length_(string_length), sum_jx2_(sum_jx2), 
        //inner_narrow_(inner_narrow), 
        //number_roots_(number_roots_total),
        lambda_(0.), epsilon_(string_length), delta_(string_length), 
        hold_(string_length, false),
        iterations_(0), 
        initial_deviation_(1e-10), 
        run_max_(50), 
        damping_delta_(0.5), 
        steps_no_deviation_(0)
{
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
    if (string_length_%2) 
        return vector<double> (1, lambda_);
    return vector<double> ();
}

    
vector< complex<double> > String::getRoots() const
{
    vector< complex<double> > result (string_length_*2);
    for (int a=0; a<string_length_; ++a) {
        // run through string, this covers conjugates.
        result[a] = lambda_ + epsilon_[a] + I*(delta_[a] + 0.5*(string_length_-1) - a);
        // opposite string
        result[string_length_*2-a-1] = -result[a];
    }
    return result;
}

int String::size() const
{
    return string_length_*2;
}



/** calculate a new value for the deviance (delta) and aberration (epsilon) **/
bool String::stepDeviation (
                const vector<SymRoots*>& all, const int alpha, 
                vector<double>& new_delta, vector<double>& new_epsilon)
{
    // outputs
	new_delta.resize(string_length_);
	new_epsilon.resize(string_length_);
	
	// real rapidity, no deviations.
	if (1==string_length_) 		{
		new_delta[0] = new_epsilon[0] = 0;
		return true;
	}

	long double theta = 0.0;
	long double log_r = 0.0;
	
	long double del_term[string_length_/2];
	long double eps_term[string_length_/2];
	
	int number_roots_odd = 0;
	for (int a=0; a<string_length_/2; ++a) {
		
		number_roots_odd = 0;
		 
		double x = lambda_ + epsilon_[a];
	    double y = 0.5*(string_length_-1) - a + delta_[a];
	    
		// string self scattering minus dangerous terms		
		for (int b=0; b < string_length_; ++b) {
			
			double x_sum = 2.*lambda_ + epsilon_[a] + epsilon_[b];
			double y_sum = (string_length_-1.) - (a+b) + delta_[a] + delta_[b];
	    
	        double x_diff = epsilon_[a] - epsilon_[b];
	        double y_diff = (delta_[a]-delta_[b]) - (a-b);
	    
			if (b!=a) {    
			    if (a-b!=1)  theta -=  xi( x_diff, 1.+y_diff );
			    if (a-b!=-1) theta -=  xi( x_diff, 1.-y_diff );

			    if (a-b!=1)  log_r +=  log( sq(x_diff) + sq(1.+y_diff) );
			    if (a-b!=-1) log_r += -log( sq(x_diff) + sq(1.-y_diff) );

		        // opposite string, excluding opposite root absorbed in kinetic term
			    theta -=  xi( x_sum, 1.-y_sum );
			    theta -=  xi( x_sum, 1.+y_sum );

			    log_r +=  log( sq(x_sum) + sq(1.+y_sum) );
			    log_r += -log( sq(x_sum) + sq(1.-y_sum) );
		    }

		}		
	    
		// scattering with others
		for (int beta=0; beta<all.size(); ++beta) {
		    // self scattering covered above
		    if (beta==alpha) continue;
		
		
            //SymRoots* roots_beta = static_cast<SymRoots*>(all[beta]);		
            vector<double> other_x = all[beta]->getRealPairs();
            vector<double> other_y = all[beta]->getImagPairs();
            vector<complex<double> > other_z = all[beta]->getComplexQuartets();
            
            if (all[beta]->hasOrigin()) {
                // there can only be one origin root
                //if (number_roots_odd) return false;
                ++number_roots_odd;
                
                theta -=  xi(x, 1.+y) 
                        + xi(x, 1.-y);
                log_r +=  log( sq(x) + sq(1.+y) )
						- log( sq(x) + sq(1.-y) );
            }
            
            for (int k=0; k<other_x.size(); ++k) {
                double x_k = other_x[k];
                
                // +x
			    theta -=  xi( x-x_k, 1.+y )
						+ xi( x-x_k, 1.-y );
				log_r +=  log( sq(x-x_k) + sq(1.+y) )
						- log( sq(x-x_k) + sq(1.-y) );
	            // -x
			    theta -=  xi( x+x_k, 1.+y )
						+ xi( x+x_k, 1.-y );
				log_r +=  log( sq(x+x_k) + sq(1.+y) )
						- log( sq(x+x_k) + sq(1.-y) );
	        }
            
            for (int k=0; k<other_y.size(); ++k) {
                double y_k = other_y[k];
                
                // +iy
			    theta -=  xi( x, 1.0 + (y-y_k) )
						+ xi( x, 1.0 - (y-y_k) );
				log_r +=  log( sq(x) + sq(1.0 + (y-y_k)) )
						- log( sq(x) + sq(1.0 - (y-y_k)) );
	            // -iy
			    theta -=  xi( x, 1.0 + (y+y_k) )
						+ xi( x, 1.0 - (y+y_k) );
				log_r +=  log( sq(x) + sq(1.0 + (y+y_k)) )
						- log( sq(x) + sq(1.0 - (y+y_k)) );
	        }
            
            for (int k=0; k<other_z.size(); ++k) {
                double x_k = real(other_z[k]);
                double y_k = imag(other_z[k]);
                // x+iy
			    theta -=  xi( x-x_k, 1.0 + (y-y_k) ) + xi( x-x_k, 1.0 - (y-y_k) );
				log_r +=  log( sq(x-x_k) + sq(1.0 + (y-y_k)) )- log( sq(x-x_k) + sq(1.0 - (y-y_k)) );
	            // -x-iy
			    theta -=  xi( x+x_k, 1.0 + (y+y_k) ) + xi( x+x_k, 1.0 - (y+y_k) );
				log_r +=  log( sq(x+x_k) + sq(1.0 + (y+y_k)) ) - log( sq(x+x_k) + sq(1.0 - (y+y_k)) );
                // x-iy
			    theta -=  xi( x-x_k, 1.0 + (y+y_k) ) + xi( x-x_k, 1.0 - (y+y_k) );
				log_r +=  log( sq(x-x_k) + sq(1.0 + (y+y_k)) ) - log( sq(x-x_k) + sq(1.0 - (y+y_k)) );
	            // -x+iy
			    theta -=  xi( x+x_k, 1.0 + (y-y_k) ) + xi( x+x_k, 1.0 - (y-y_k) );
				log_r +=  log( sq(x+x_k) + sq(1.0 + (y-y_k)) ) - log( sq(x+x_k) + sq(1.0 - (y-y_k)) );
	        
	        }
	    }
        
        
		// contribution M mod 2 from sum of conjugate quantum numbers
		theta -= PI*number_roots_odd; //PI*number_roots_		
		
		// kinetic theta, absorbing opposite root
		theta += (big_n_-1.) * ( xi(2.*x, 1.+2.*y) + xi(2.*x, 1.-2.*y) );
		log_r += (big_n_-1.) * ( log(sq(x) + sq(y-0.5)) - log(sq(x) + sq(y+0.5)) );

		// theta is now sum b=1..a theta[b]
		// log_r is now sum b=1..a log_r[b]
		del_term[a] = -cos(theta)*exp(0.5*log_r);
		eps_term[a] = sin(theta)*exp(0.5*log_r);
	}	
		
    
	if (string_length_%2) {
	    int a_mid = (string_length_-1)/2;
	    
		// odd n_j
		for (int a=0; a < a_mid; ++a) {
			for (int b=a; b < string_length_/2; ++b) {
				new_delta[a]   += del_term[b];
				new_epsilon[a] += eps_term[b];
			}
		}
		
		// string center has no deviations.
		new_delta[a_mid] = 0.0;
		new_epsilon[a_mid] = 0.0;
		
		// set the conjugates
		for (int a=a_mid+1; a<string_length_;++a) {
			new_delta[a] = - new_delta[string_length_-a-1];
			new_epsilon[a] = new_epsilon[string_length_-a-1];
		}
	}
	else {
		// even n_j. 
		// this is more involved, see notes
	    int a_mid = string_length_/2 - 1;
	    	
		// first, the regular terms.
		for (int a=0; a<a_mid; ++a) {
			for (int b=a; b<a_mid; ++b) {
				new_delta[a] += del_term[b];
				new_epsilon[a] += eps_term[b];
			}
		}

        // apply the inner-sign term.
        bool inner_narrow = innerPairIsNarrow(sum_jx2_, string_length_, number_roots_odd);
		for (int a=0; a<a_mid+1; ++a) {
			// log_r should be sum b=1...n/2 log_r[b] now
			new_delta[a] += 0.5*exp(0.5*log_r) * (inner_narrow?-1.:1.); 
		}
		
		// inner pair has no aberration
		new_epsilon[a_mid] = 0;
		
		// set the conjugates
		for (int a=a_mid+1; a<string_length_; ++a) {
			new_delta[a] = - new_delta[string_length_-a-1];
			new_epsilon[a] = new_epsilon[string_length_-a-1];
		}
	}
	
	return true;
}



/** calculate a new value for the string centre (lambda) **/
double String::stepRapidity (const vector<SymRoots*>& all, const int alpha)
{
	long double scatter = 0.0;		
	for (int a=0; a < string_length_; ++a) {

        double x = lambda_ + epsilon_[a];
        double y = 0.5*(string_length_-1) - a + delta_[a];
            
        // scattering with opposite strings
	    for (int b=0; b < string_length_; ++b) {
            // direct opposite - absorb into kinetic term
            if (b+a==string_length_-1) continue;
 
            double x_b = - lambda_ - epsilon_[b];
	        double y_b = 0.5*(string_length_-1) - b + delta_[b];
 			scatter += xi( x-x_b, 1.-(y-y_b) )  +  xi( x-x_b, 1.+(y-y_b) );
		}

	
		// scattering with others
		for (int beta=0; beta<all.size(); ++beta) {

		    if (beta==alpha) continue;
            
            //SymRoots* roots_beta = static_cast<SymRoots*>(all[beta]);		
            vector<double> other_x = all[beta]->getRealPairs();
            vector<double> other_y = all[beta]->getImagPairs();
            vector<complex<double> > other_z = all[beta]->getComplexQuartets();

            if (all[beta]->hasOrigin()) {
				scatter += xi( x, 1.-y ) + xi( x, 1.+y );
		    }    		   

    	    for (int k=0; k < other_x.size();++k) {
    	        double x_k = other_x[k];
				scatter += xi( x-x_k, 1.-y ) + xi( x-x_k, 1.+y );
				scatter += xi( x+x_k, 1.-y ) + xi( x+x_k, 1.+y );
			}
                	
    	    for (int k=0; k < other_y.size();++k) {
    	        double y_k = other_y[k];
				scatter += xi( x, 1.-(y-y_k) )  +  xi( x, 1.+(y-y_k) );
				scatter += xi( x, 1.-(y+y_k) )  +  xi( x, 1.+(y+y_k) );
			}
    	    
    	    for (int k=0; k < other_z.size();++k) {
    	        double x_k = real(other_z[k]);
    	        double y_k = imag(other_z[k]);
				scatter += xi( x-x_k, 1.-(y-y_k) )  +  xi( x-x_k, 1.+(y-y_k) );
				scatter += xi( x+x_k, 1.-(y-y_k) )  +  xi( x+x_k, 1.+(y-y_k) );
				scatter += xi( x-x_k, 1.-(y+y_k) )  +  xi( x-x_k, 1.+(y+y_k) );
				scatter += xi( x+x_k, 1.-(y+y_k) )  +  xi( x+x_k, 1.+(y+y_k) );
			}

		}
	}
	
	// 2\pi \sum_{a=1}^n J_a + \theta\sub{other}^a = \sum_{a=1}^n \theta\sub{kin}^a 
    //long double phase = 0.5*(PI*sum_jx2_ + scatter) / N_;
    
    // absorbed opposite root in N-1 for faster convergence
    long double phase = 0.5*(PI*sum_jx2_ + scatter) / (big_n_-1.);

	double new_lambda = 0.;
	if (string_length_%2) {
		// odd length
		long double sum_kin = 0.0;
		for (int a=0; a< (string_length_-1)/2; ++a) {
			sum_kin += 	  xi( lambda_ + epsilon_[a+1], 0.5*string_length_ - (a+1) + delta_[a+1]  )
						+ xi( lambda_ + epsilon_[a],  -0.5*string_length_ + (a+1) - delta_[a]  );
		}
		new_lambda = -epsilon_[0] + ( 0.5*string_length_ + delta_[0] )*tan( phase-sum_kin ) ;
	}
	else {
		// even length
		long double sum_kin = 0.0;
		for (int a=0; a< string_length_/2-1; ++a) {
			sum_kin += 	  xi( lambda_ + epsilon_[a+1], 0.5*string_length_ - (a+1) + delta_[a+1] )
						+ xi( lambda_ + epsilon_[a],  -0.5*string_length_ + (a+1) - delta_[a] );
		}
		sum_kin += xi( lambda_, -delta_[string_length_/2-1]); 
		new_lambda = -epsilon_[0] + ( 0.5*string_length_ + delta_[0] )*tan( phase-sum_kin ) ;
    

/*
		double tan_theta =  tan(  phase - sum_kin  ) ;
		double eps_1 = epsilon_[0];
		double del_1 = delta_[0];
		double del_n2 = delta_[string_length_/2-1];
		
		double b = (eps_1*tan_theta + del_1 - del_n2 + 0.5*string_length_);
		double c = del_n2*(tan_theta*(0.5*string_length_+del_1) - eps_1);
		// keep the sign as it was
		//int sign_sqrt = isgn(2.0*lambda_*tan_theta + b);
		new_lambda = (-b + sign_sqrt*sqrt( sq(b) - 4.0*tan_theta*c ))/(2.0*tan_theta);
*/

	}

	return new_lambda;
}



bool String::iterate(const vector<SymRoots*>& all, const int alpha, const double last_convergence, double& convergence)
{
	
	//double square_diff = 0.0;	
	//double one_over_n = 1.0/big_n_;
	long double new_rap = 0.0;
	double current_rapidity = lambda_;	

	for (int run = 0; run < run_max_; ++run) {
		long double last_rap = new_rap;
		new_rap = stepRapidity(all, alpha);

		if (finite(new_rap)) break;
		// if it doesn't work, slowly raise the offending rapidity and see if we can get something reasonable.
		lambda_ *= 1.1;
	
		// growing over the machine maximum or otherwise non-finite? give up.
		if (!finite(lambda_)) break; 
	}
	lambda_ = current_rapidity; 

	// if not on hold and not finite, there's nothing we can do.
	if (!finite(new_rap)) return false;
	    //throw Exception (here, exc_Runaway);

	// calculate convergence measure
	convergence += sq(lambda_-new_rap);
    
    vector<double> new_delta(delta_.size(), 1e-20);
	vector<double> new_epsilon(epsilon_.size(), 1e-20);
	 
    if (iterations_ > steps_no_deviation_) {
	    // find deviations & calculate their convergence. if not found, give up.
	    if (!stepDeviation(all, alpha, new_delta, new_epsilon)) return false;

        // apply damping to delta
        for (int a=0; a<string_length_; ++a) 
	        new_delta[a] = (1.0-damping_delta_)*new_delta[a] + damping_delta_*delta_[a];
    }
    
    // enforce holds
    for (int a=0; a<string_length_;++a) {
	    if (hold_[a]) {
		    new_delta[a] = delta_[a];
		    new_epsilon[a] = epsilon_[a];
	    }
    }

    // calculate convergence
	for (int a=0; a < string_length_; ++a)
		convergence += sq(new_delta[a] - delta_[a]) + sq(new_epsilon[a] - epsilon_[a]);

	// set the new values for the rapidities
	lambda_ = new_rap;
	delta_ = new_delta;
	epsilon_ = new_epsilon;

    ++iterations_;
    
	// report success	
	return true;
}

	
bool String::initiate(const vector<SymRoots*>& all, const int alpha)
{
        


    if (sum_jx2_ == 0)
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
    
    for(int a=0; a< string_length_/2; ++a) {
        epsilon_[a] = initial_deviation_;
        delta_[a] = initial_deviation_;
    }
    for (int a=string_length_-string_length_/2; a<string_length_; ++a) {
        epsilon_[a] = initial_deviation_;
        delta_[a] = -initial_deviation_;
    }
		
    // for odd strings, set central deviations.
    // for even strings, ensure the deviation of the inner pair has the correct sign
    if (string_length_%2) {
        delta_[string_length_/2] = 0.;
        epsilon_[string_length_/2] = initial_deviation_;
    }
    else {
        bool number_roots_odd = false;
        for (int beta=0; beta<all.size();++beta) { 
            //SymRoots* roots_beta = static_cast<SymRoots*>(all[beta]); 
            if (all[beta]->hasOrigin()) {
                number_roots_odd=true;
                break;
            } 
        }
        if (innerPairIsNarrow(sum_jx2_, string_length_, number_roots_odd)) {
            delta_[string_length_/2-1] = -initial_deviation_;
            delta_[string_length_/2] = initial_deviation_;       
        } 
    }
    
    return true;
}







/* CentralString */


CentralString::CentralString(const int big_n, const int string_length)
      : big_n_(big_n), string_length_(string_length), sum_jx2_(0), 
        number_roots_odd_(false),
        delta_(string_length), 
        hold_(string_length, false), knockout_(0),
        iterations_(0)
{
    // use default policies, see header file
    setPolicy();
}

void CentralString::setPolicy(const double initial_deviation, const double damping_delta)
{
    initial_deviation_ = initial_deviation;
    damping_delta_ = damping_delta;
}    

bool CentralString::knockOut(const int number_pairs)
{
    if (0==number_pairs) knockout_ = 0; 
    else knockout_ = 2*number_pairs - (string_length_%2);
    return (knockout_ <= string_length_);
} 


vector<double> CentralString::getImagPairs() const
{
    // run through half string, excluding odd strings' centres 
    const int n_pairs = (string_length_-knockout_)/2;
    vector<double> result(n_pairs);
    for (int a=0; a<n_pairs; ++a) {
        result[a] = delta_[a] + 0.5*(string_length_-1) - a;
    }
    return result;
}

vector<double> CentralString::getImagPairsDelta() const
{
    // run through half string, excluding odd strings' centres
    const int n_pairs = (string_length_-knockout_)/2 ;
    vector<double> result(n_pairs);
    for (int a=0; a<n_pairs; ++a) {
        result[a] = delta_[a];
    }
    return result;
}

vector<int> CentralString::getImagPairsPos() const
{
    // run through half string, excluding odd strings' centres
    const int n_pairs = (string_length_-knockout_)/2 ;
    vector<int> result(n_pairs);
    for (int a=0; a<n_pairs; ++a) {
        result[a] = (string_length_-1) - 2*a;
    }
    return result;
}

bool CentralString::hasOrigin() const
{
    return (!knockout_ && string_length_%2);
}

    

vector< complex<double> > CentralString::getRoots() const
{
    const int n_roots = string_length_ - knockout_;
    vector< complex<double> > result (n_roots);
    for (int a=0; a< (string_length_-knockout_)/2;  ++a) {
        result[a] = I*(delta_[a] + 0.5*(string_length_-1) - a);
        result[n_roots-1-a] = conj(result[a]);
    }
    if (!knockout_ && string_length_%2) 
        result[string_length_/2] = 0.; 
    return result;
}


int CentralString::size() const
{
    return string_length_ - knockout_;
}



bool CentralString::iterate(const vector<SymRoots*>& all, const int alpha, const double last_convergence,double& convergence)
{
	vector<double> step_delta;
	// find deviations & calculate their convergence. if not found, give up.
    if (!stepCentralDeviation(all, alpha, step_delta)) return false;

    // secant method step
    vector<double> new_delta (string_length_);
    const int a_last = knockout_? (string_length_-knockout_)/2 : (string_length_-1)/2;  //exclude central pair even if no knockout 
    const int a_mid = (string_length_-1)/2; 
    
    for (int a=0; a<a_last; ++a)  {
        new_delta[a] = secantStep(delta_last_[a], delta_[a], step_last_[a], step_delta[a], convergence);
        //new_delta[a] = step_delta[a];
        new_delta[string_length_-1-a] = - new_delta[a];
    }
    for (int a=a_last; a<=a_mid; ++a)  {
        new_delta[a] = step_delta[a];
        new_delta[string_length_-1-a] = - new_delta[a];
    }    

/*        
    // damping
    for (int a=0; a<string_length_; ++a) 
        new_delta[a] = (1.0-damping_delta_)*new_delta[a] + damping_delta_*delta_[a];
*/
    
    // enforce holds
    for (int a=0; a<string_length_; ++a) 
	    if (hold_[a]) new_delta[a] = delta_[a];
    
	
    // calculate convergence
	//for (int a=0; a < string_length_; ++a)
	//	convergence += sq(new_delta[a] - delta_[a]);

    delta_last_ = delta_;
    step_last_ = step_delta;
	delta_ = new_delta;
    
    ++iterations_;
	// report success	
	return true;
}


void CentralString::setInitialValues(const double initial_deviation) //const bool number_roots_odd, 
{
    for(int a=0; a<string_length_/2; ++a) {
        delta_[a] = initial_deviation * double(string_length_/2-a);
        delta_[string_length_-1-a] = -delta_[a];
    }
    // sets central root/inner pair to zero
    delta_[(string_length_-1)/2] = 0.;
    delta_[string_length_/2] = 0.;
}    

	
bool CentralString::initiate(const vector<SymRoots*>& all, const int alpha)
{
    number_roots_odd_ = false;
	// run over all complexes, including ourselves
	for (int beta=0; beta<all.size(); ++beta) {
	    if (all[beta]->hasOrigin()) { 
	        // there can (should) be only one origin root
	        // and number or roots is odd if and only if it's there
	        number_roots_odd_ = true;
	        break;
	    }
	}    
	
	// set initial values to prepare for secant mehod    
    setInitialValues(initial_deviation_);
    if (!stepCentralDeviation(all, alpha, step_last_)) return false;
    delta_last_ = delta_;
    setInitialValues(2.*initial_deviation_);
    
    return true;
}



void CentralString::scatterOthers(const vector<SymRoots*>& all, int& theta, long double& log_rsq, const int alpha, const int a) const
{

    const long double y = 0.5*(string_length_-1) - a + delta_[a];
    const long double del = delta_[a];
    const int pos = (string_length_-1) - 2*a;

	// scattering with others
	for (int beta=0; beta<all.size(); ++beta) {
	    // self scattering covered above
	    if (beta==alpha) continue;

        const vector<double> other_x = all[beta]->getRealPairs();
        const vector<double> other_del = all[beta]->getImagPairsDelta();
        const vector<int> other_pos = all[beta]->getImagPairsPos();
        const vector<complex<double> > other_z = all[beta]->getComplexQuartets();
        
        if (all[beta]->hasOrigin()) {
            // there should only be one origin root
            //++number_roots_odd;
            
            theta -=  (abs(y)>1.) ? sgn(y) : 0;
            log_rsq +=  4.L*atih(y); //log( sq((1.L+y) / (1.L-y)) );
            
            // 4.L* --> 2* because we're looking at log r^2 rather than log_r, another 2* because we've added eq and conj(eq)
        }

        for (int k=0; k<other_x.size(); ++k) {
            const long double x_k = other_x[k];
            // +x // -x
			log_rsq +=  2.L*log( (sq(x_k) + sq(1.L+y)) / (sq(x_k) + sq(1.L-y)) );
        }

        for (int k=0; k<other_del.size(); ++k) {
            // split position and delta for precision in case positions are equal
            const long double del_dif = del-other_del[k];
            const long double del_sum = del+other_del[k];
            const long double pos_dif = 0.5L*(pos-other_pos[k]);
            const long double pos_sum = 0.5L*(pos+other_pos[k]);
            
            // +iy
		    theta -=  (abs(pos_dif+del_dif)>1.) ? sgn(pos_dif+del_dif) : 0;
			log_rsq +=  4.L*atih(pos_dif+del_dif);
            // -iy
		    theta -=  (abs(pos_sum+del_sum)>1.) ? sgn(pos_sum+del_sum) : 0;
			log_rsq +=  4.L*atih(pos_sum+del_sum);
        }

        for (int k=0; k<other_z.size(); ++k) {
            const long double x_k = real(other_z[k]);
            const long double y_k = imag(other_z[k]);
            // x+iy  -x+iy
            log_rsq +=  2.L*(log( sq(x_k) + sq(1.L + (y-y_k)) ) - log( sq(x_k) + sq(1.L - (y-y_k)) ));
            // -x-iy  x-iy
			log_rsq +=  2.L*(log( sq(x_k) + sq(1.L + (y+y_k)) ) - log( sq(x_k) + sq(1.L - (y+y_k)) ));
        }
    }
    return;       
}


/** calculate a new value for the deviance (delta) and aberration (epsilon) **/
bool CentralString::stepCentralDeviation (
                const vector<SymRoots*>& all, const int alpha, 
                vector<double>& new_delta)
{
    // real rapidity, no deviations.
	if (1==string_length_) 	{
	    new_delta = { 0. };
		return true;
	}

	long double sum_log_rsq = 0.L;
	int sum_theta = 0; // number of PI terms
	vector<long double> del_term ((string_length_+1)/2, 0.);
	long double log_rsq_other = 0.L;


    const int a_last = (string_length_-knockout_-1)/2;    
	
	for (int a=0; a <= a_last; ++a) {
	    
        {
            long double log_rsq = 0.L;
            // scattering with others
            scatterOthers(all, sum_theta, log_rsq, alpha, a);
            // keep for use in sum equation
            log_rsq_other += log_rsq;
            sum_log_rsq += log_rsq;
        }
        if (a==a_last) break;
        
        // kinetic term
		const long double y = 0.5*(string_length_-1) - a + delta_[a];
		// kinetic theta 
		sum_theta += (abs(2.*y)>1.) ? big_n_*sgn(2.*y) : 0;
		// 4.L*atih(2.L*y) == log( sq((y+0.5L) / (y-0.5L)) ) 
		sum_log_rsq -= (long double)big_n_ * 4.L*atih(2.L*y); 
		
        // contribution from 2*quantum number
		sum_theta -= 1-number_roots_odd_;


		// string self scattering minus dangerous terms		
		for (int b=0; b < string_length_; ++b) {
		    // skip knocked out roots
		    if ( b>=(string_length_-knockout_)/2 && b<(string_length_+knockout_)/2 ) continue;

			// no need to split pos and delta here as we're excluding dangerous terms
	        const long double y_diff = (delta_[a]-delta_[b]) - (a-b);
			if (b!=a) {    
                if (a-b!=+1 && a-b!=-1)  sum_theta -= (abs(y_diff)>1.) ? sgn(y_diff) : 0;
			    if (a-b!=+1)  sum_log_rsq +=  log(sq( 1.L+y_diff ));
			    if (a-b!=-1)  sum_log_rsq -=  log(sq( 1.L-y_diff ));
		    }
		}	
            
	    const long double sign = (sum_theta%2) ? 1: -1;
		del_term[a] = sign*exp(0.5L*sum_log_rsq);  
	}	


    new_delta.clear();
	new_delta.resize(string_length_, 0.);

    // I think formally correct, but
    // does not converge
    if (false) {
        if(knockout_) {
                
            long double new_delta_last = 0.;
            long double sum_kin = 0.0;
		
            // log(sq(2.L*y+1.L))- log(sq(2.L*y-1.L));
            // other kinetic terms, paired to form logsqs close to zero (factors close to one) 
	        for (int a=0; a< a_last; ++a) {
		        sum_kin += log(sq(   ( (string_length_-2*(a+1)) + 2.L*delta_[a+1] )
		                            /(-(string_length_-2*(a+1)) - 2.L*delta_[a]   ) ));
	        }

    long double y0 = 0.5*(string_length_-1) - 0 + delta_[0];
    long double y1 = 0.5*(string_length_-1) - 1 + delta_[1];
    long double pos1 = (0.5*(string_length_-1) - 1);    
    cerr<<"? "<< 0.5* (-2.*pos1-1. + ((2.*y0+1.))/exp(-0.5*sum_kin+ 0.5*log_rsq_other/big_n_))<<endl;
    cerr<< 0.5*(-2.*pos1-1. + ((2.*y0+1.))/(  (2.*y0-1.)/(2.*y1+1.)*exp(+0.5*log_rsq_other/big_n_)))<<endl;
    cerr<<"? "<< 0.5* ( -string_length_+2.) + 0.5*( string_length_ + 2.*delta_[0] ) / exp(-0.5*sum_kin+ 0.5*log_rsq_other/big_n_)<<endl;
    cerr<<"?? "<<-0.5*(knockout_+2) + 0.5*(string_length_ +2.L*delta_[0]) / exp(0.5L*(-sum_kin + log_rsq_other/big_n_))<<endl;
		
            // missing terms: 
            //  + logsq( string_length_ +2.L*delta[0] ) 
            //  - logsq ( - knockout_ - 2.L*delta[a_last])  [if l even or knockout!=0, so in this block always]  
	        new_delta_last = -0.5*(knockout_+2) + 0.5*(string_length_ +2.L*delta_[0]) / exp(0.5L*(-sum_kin + log_rsq_other/big_n_));
cerr<<new_delta_last<<endl;	    	
	        new_delta[a_last] = new_delta_last;
	    }
	    
	    // add differences
	    for (int a=a_last-1; a >= 0; --a)   
		    new_delta[a] += new_delta[a+1] + del_term[a];
	}    
	
	
	// this version works for a number of cases (where a_last=0)
	// formally incorrect in other cases, and doesn't converge
    if (false) {
        if (knockout_) {
            long double new_delta_last = 0.;
            long double sum_kin = 0.0;
		
            // log(sq(2.L*y+1.L))- log(sq(2.L*y-1.L));
            // other kinetic terms, paired to form logsqs close to zero (factors close to one) 
	        for (int a=0; a< a_last; ++a) {
		        // kinetic term
		        const long double y = 0.5*(string_length_-1) - a + delta_[a];
		        // 4.L*atih(2.L*y) == log( sq((y+0.5L) / (y-0.5L)) ) 
		        sum_kin += 4.L*atih(2.L*y); 
            }

            new_delta_last =0.5/tanh( (-sum_kin + log_rsq_other) / (4.*(big_n_-1.)) ) - 0.5*(string_length_-1) + a_last;
	        new_delta[a_last] = new_delta_last;
	    }
	    
	    // add differences
	    for (int a=a_last-1; a >= 0; --a)   
		    new_delta[a] += new_delta[a+1] + del_term[a];
	
	}    
    
    // also works for a_last=0 only
    // formally incorrect in other cases, and doesn't converge
    
    if (false) {
        if (knockout_) {
            
            long double sum_kin = 0.0;
		
            // log(sq(2.L*y+1.L))- log(sq(2.L*y-1.L));
            // other kinetic terms, paired to form logsqs close to zero (factors close to one) 
	        for (int a=1; a<= a_last; ++a) {
		        // kinetic term
		        const long double y = 0.5*(string_length_-1) - a + delta_[a];
		        // 4.L*atih(2.L*y) == log( sq((y+0.5L) / (y-0.5L)) ) 
		        sum_kin += 4.L*atih(2.L*y); 
            }

	        new_delta[0] = 0.5/tanh( (-sum_kin + log_rsq_other) / (4.*(big_n_-1.)) ) - 0.5*(string_length_-1) + 0;
	    }
	    
	    // add differences
	    for (int a=1; a <= a_last; ++a)   
		    new_delta[a] += new_delta[a-1] - del_term[a-1];
	}    


    // i think thos one is at least formally correct, but again no convergence except if a_last=0
    if (true) {
        if (knockout_) {
            
            long double sum_kin = 0.0;
		
            // log(sq(2.L*y+1.L))- log(sq(2.L*y-1.L));
            // other kinetic terms, paired to form logsqs close to zero (factors close to one) 
	        
	        { int a = 0;
		        // kinetic term
		        //const long double y = 0.5*(string_length_-1) - a + delta_[a];
		        // 4.L*atih(2.L*y) == log( sq((y+0.5L) / (y-0.5L)) ) 
		        //sum_kin += 4.L*atih(2.L*y); 
            
                log_rsq_other = 0;
                int theta = 0;
                scatterOthers(all, theta, log_rsq_other, alpha, a);

		        // string self scattering 		
		        for (int b=0; b < string_length_; ++b) {
			        // skip knocked out roots
		            if ( b>=(string_length_-knockout_)/2 && b<(string_length_+knockout_)/2 ) continue;

	                const double y_diff = (delta_[a]-delta_[b]) - (a-b);
	            
			        if (b!=a && b != string_length_-1-a) {    
			            // split position and delta differences for precision in the case where (1-a+b)==0 etc.
            	        log_rsq_other += log(sq( ((double)(1-a+b) + (delta_[a]-delta_[b])) / ((double)(1+a-b) - (delta_[a]-delta_[b])) ));
                    }   
                }		    
            }

	        new_delta[0] = 0.5/tanh( (-sum_kin + log_rsq_other) / (4.*(big_n_-1.)) ) - 0.5*(string_length_-1) + a_last;
	    }
	    
	    // add differences
	    for (int a=1; a <= a_last; ++a)   
		    new_delta[a] += new_delta[a-1] - del_term[a-1];
	}    

	
	// l=5: 0 1 2 3 4, a_mid=2
    // l=4: 0 1 2 3, a_mid=1
    const int a_mid = (string_length_-1)/2;
    // set the conjugates
	// note that inner pair/string centre are necessarily zero - 
    // have been zeroed on initialisation 
    // new_delta[a_mid] = 0.;
	for (int a=a_mid+1; a<string_length_;++a) 
		new_delta[a] = - new_delta[string_length_-a-1];
	
	return true;
}





vector<complex<double> > CentralString::getCleanerJ (
                            const vector<SymRoots*>& all, const int alpha) const 
{
    // outputs
	vector<complex<double> > result(string_length_);
	
	for (int a=0; a<string_length_; ++a) {
		long double log_rsq = 0.0;
	    int theta = 0;
	
		const double y = 0.5*(string_length_-1) - a + delta_[a];
        const double del = delta_[a];
	    const int pos = (string_length_-1) - 2*a;
	    
		// string self scattering 		
		for (int b=0; b < string_length_; ++b) {
			
	        const double y_diff = (delta_[a]-delta_[b]) - (a-b);
	    
			if (b!=a) {    
			    theta -=  (abs(y_diff)>1.) ? sgn(y_diff) : 0;
			    // split position and delta differences for precision in the case where (1-a+b)==0 etc.
    	        log_rsq += log(sq( ((double)(1-a+b) + (delta_[a]-delta_[b])) / ((double)(1+a-b) - (delta_[a]-delta_[b])) ));
            }   
        }		
		
        scatterOthers(all, theta, log_rsq, alpha, a);
            
        
		// kinetic theta 
		theta += (abs(2.*y)>1.) ? big_n_*sgn(2.*y) : 0;
		log_rsq += (big_n_) * log( sq(y-0.5) / sq(y+0.5));

        // imaginary part should be small, is an error term.
        result[a] = complex<double> (0.5*theta, 0.5*log_rsq);      
    }	
		
    return result;
}

// error is given as average of square error per quantum number
vector<int> CentralString::getCleanJx2 (
                        const vector<SymRoots*>& all, const int alpha, double& sq_error) const
{
    const vector< complex<double> > dirty_j = getCleanerJ(all, alpha);
    vector<int> result(dirty_j.size(), 0);
    sq_error=0.;
    for (int i=0; i<result.size(); ++i) {
        result[i] = round(2.*real(dirty_j[i]));
        sq_error += norm(1.*result[i] - 2.*dirty_j[i])/(double)dirty_j.size();           
    }      
    return result;
}
