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


CentralString::CentralString(const int big_n, const int string_length)// /*const int number_roots_total, const bool inner_narrow*/
      : big_n_(big_n), string_length_(string_length), sum_jx2_(0), 
         delta_(string_length), 
        hold_(string_length, false),
        iterations_(0), 
        initial_deviation_(1e-10), 
        damping_delta_(0.5) 
{
}


 bool CentralString::knockOut(const int number_pairs) {};


vector< double > CentralString::getImagPairs() const
{
    vector< double > result(string_length_/2);
    // run through half string, excluding odd strings' centres 
    for (int a=0; a<(string_length_)/2; ++a) {
        result[a] =  (delta_[a] + 0.5*(string_length_-1) - a);
    }
    return result;
}

vector< double > CentralString::getImagPairsDelta() const
{
    vector< double > result(string_length_/2);
    // run through half string, excluding odd strings' centres 
    for (int a=0; a<string_length_/2; ++a) {
        result[a] = delta_[a];
    }
    return result;
}

vector< int > CentralString::getImagPairsPos() const
{
    vector< int > result(string_length_/2);
    // run through half string, excluding odd strings' centres 
    for (int a=0; a<string_length_/2; ++a) {
        result[a] = (string_length_-1) - 2*a;
    }
    return result;
}



bool CentralString::hasOrigin() const
{
    return (string_length_%2);
}
    

vector< complex<double> > CentralString::getRoots() const
{
    vector< complex<double> > result (string_length_);
    for (int a=0; a<string_length_; ++a) {
        // run through string, this covers conjugates.
        result[a] = I*(delta_[a] + 0.5*(string_length_-1) - a);
        // no opposite string
    }
    return result;
}

int CentralString::size() const
{
    return string_length_;
}




bool CentralString::iterate(const vector<SymRoots*>& all, const int alpha, const double last_convergence,double& convergence)
{
	//double one_over_n = 1.0/big_n_;
	
	vector<double> new_delta;
    // find deviations & calculate their convergence. if not found, give up.
    if (!stepCentralDeviation(all, alpha, new_delta)) return false;
    // apply damping to delta
    for (int a=0; a<string_length_; ++a) 
        new_delta[a] = (1.0-damping_delta_)*new_delta[a] + damping_delta_*delta_[a];

    
    // enforce holds
    for (int a=0; a<string_length_;++a) {
	    if (hold_[a]) {
		    new_delta[a] = delta_[a];
		}
    }
    
    // enforce pm 0.5i
    //if (!(string_length_%2)) {
       //new_delta[string_length_/2-1] = 0.;
       //new_delta[string_length_/2] = 0.;
    //}
	
    // calculate convergence
	for (int a=0; a < string_length_; ++a)
		convergence += sq(new_delta[a] - delta_[a]);

	delta_ = new_delta;
    ++iterations_;
	// report success	
	return true;
}

	
bool CentralString::initiate(const vector<SymRoots*>& all, const int alpha)
{
    
    for(int a=0; a< string_length_/2; ++a) {
        delta_[a] = initial_deviation_ * double(string_length_/2-a);
    }
    for (int a=string_length_-string_length_/2; a<string_length_; ++a) {
        delta_[a] = -initial_deviation_ * double(a);
    }
		
    // for odd strings, set central deviations.
    // for even strings, ensure the deviation of the inner pair has the correct sign
    if (string_length_%2) {
        delta_[string_length_/2] = 0.;
    }
    else {
        bool number_roots_odd = 0;
        for (int beta=0; beta<all.size();++beta) { 
            //SymRoots* roots_beta = static_cast<SymRoots*>(all[beta]); 
            if (all[beta]->hasOrigin()) {
                number_roots_odd=true;
                break;
            } 
        }
        if (innerPairIsNarrow(0, string_length_, number_roots_odd)) {
            
            delta_[string_length_/2-1] = -initial_deviation_;
            delta_[string_length_/2] = initial_deviation_;       
        } 
    }
    
    return true;
}




/** calculate a new value for the deviance (delta) and aberration (epsilon) **/
bool CentralString::stepCentralDeviation (
                const vector<SymRoots*>& all, const int alpha, 
                vector<double>& new_delta)
{
    // outputs
    new_delta.clear();
	new_delta.resize(string_length_);
	
	// real rapidity, no deviations.
	if (1==string_length_) 		{
		new_delta[0] =  0;
		return true;
	}

	long double log_rsq = 0.L;
	int theta = 0;
	
	vector<long double> del_term (string_length_/2,0.);
	
	int number_roots_odd = 0;
	
	for (int a=0; a<string_length_/2; ++a) {
		number_roots_odd = (string_length_%2);
	
        long double y = 0.5*(string_length_-1) - a + delta_[a];
	    long double del = delta_[a];
	    int pos = (string_length_-1) - 2*a;
	    
	    
		// string self scattering minus dangerous terms		
		for (int b=0; b < string_length_; ++b) {
		
			// no need to split pos and delta here as we're excluding dangerous terms
	        long double y_diff = (delta_[a]-delta_[b]) - (a-b);
	        
			if (b!=a) {    
                if (a-b!=1 && a-b!=-1)  theta -=  (abs(y_diff)>1.) ? sgn(y_diff) : 0;
			    
			    if (a-b!=+1)  log_rsq +=  log(sq( 1.L+y_diff ));
			    if (a-b!=-1)  log_rsq -=  log(sq( 1.L-y_diff ));
		    }

		}		
        
cerr<<theta<<endl;
		
		
		// scattering with others
		for (int beta=0; beta<all.size(); ++beta) {
		    // self scattering covered above
		    if (beta==alpha) continue;
    
    		//SymRoots* roots_beta = static_cast<SymRoots*>(all[beta]); 
            vector<double> other_x = all[beta]->getRealPairs();
            vector<double> other_del = all[beta]->getImagPairsDelta();
            vector<int> other_pos = all[beta]->getImagPairsPos();
            vector<complex<double> > other_z = all[beta]->getComplexQuartets();
            
            if (all[beta]->hasOrigin()) {
                // there should only be one origin root
                ++number_roots_odd;
                
                theta -=  (abs(y)>1.) ? sgn(y) : 0;
                log_rsq +=  4.L*atih(y); //log( sq((1.L+y) / (1.L-y)) );
            }
            
            for (int k=0; k<other_x.size(); ++k) {
                long double x_k = other_x[k];
                // +x
				//log_rsq +=  log( (sq(-x_k) + sq(1.L+y)) / (sq(-x_k) + sq(1.L-y)) );
	            // -x
				//log_rsq +=  log( (sq(+x_k) + sq(1.L+y)) / (sq(+x_k) + sq(1.L-y)) );
	            log_rsq +=  2.L*log( (sq(x_k) + sq(1.L+y)) / (sq(x_k) + sq(1.L-y)) );
	        
	        }

            for (int k=0; k<other_del.size(); ++k) {
                
                // split position and delta for precision in case positions are equal
                long double del_dif = del-other_del[k];
                long double del_sum = del+other_del[k];
                
                long double pos_dif = 0.5L*(pos-other_pos[k]);
                long double pos_sum = 0.5L*(pos+other_pos[k]);
                
                // +iy
			    theta -=  (abs(pos_dif+del_dif)>1.) ? sgn(pos_dif+del_dif) : 0;
				log_rsq +=  4.L*atih(pos_dif+del_dif);
	            //log_rsq +=  log(sq( (1.L + (y-y_k)) / (1.L - (y-y_k)) ));
	            
	            // -iy
			    theta -=  (abs(pos_sum+del_sum)>1.) ? sgn(pos_sum+del_sum) : 0;
				log_rsq +=  4.L*atih(pos_sum+del_sum);
				//log_rsq +=  log(sq( (1.L + (y+y_k)) / (1.L - (y+y_k)) ));
            }
 
                   
            for (int k=0; k<other_z.size(); ++k) {
                long double x_k = real(other_z[k]);
                long double y_k = imag(other_z[k]);
                // x+iy
				//log_rsq +=  log( sq(-x_k) + sq(1.L + (y-y_k)) ) - log( sq(-x_k) + sq(1.L - (y-y_k)) );
	            // -x-iy
				//log_rsq +=  log( sq(+x_k) + sq(1.L + (y+y_k)) ) - log( sq(+x_k) + sq(1.L - (y+y_k)) );
                // x-iy
				//log_rsq +=  log( sq(-x_k) + sq(1.L + (y+y_k)) ) - log( sq(-x_k) + sq(1.L - (y+y_k)) );
	            // -x+iy
				//log_rsq +=  log( sq(+x_k) + sq(1.L + (y-y_k)) ) - log( sq(+x_k) + sq(1.L - (y-y_k)) );
	            // x-iy
				log_rsq +=  2.*(log( sq(-x_k) + sq(1.L + (y+y_k)) ) - log( sq(-x_k) + sq(1.L - (y+y_k)) ));
	            // -x+iy
				log_rsq +=  2.*(log( sq(+x_k) + sq(1.L + (y-y_k)) ) - log( sq(+x_k) + sq(1.L - (y-y_k)) ));
	        
	        }
        
cerr<<theta<<endl;

	    }
        
cerr<<theta<<endl;
        
		// kinetic theta 
		theta += (abs(2.*y)>1.) ? big_n_*sgn(2.*y) : 0;
		log_rsq -= (long double)big_n_ * 4.L*atih(2.L*y);
		//log_rsq -= (long double)(big_n) * log( sq((y+0.5L) / (y-0.5L)) ) ;
cerr<<theta<<endl;    
        // contribution from 2*quantum number
		theta -= 1-number_roots_odd;
cerr<<theta<<endl;
		
	    long double sign = (theta%2) ? 1: -1;
		del_term[a] = sign*exp(0.5L*log_rsq);  
	}	
cerr<<del_term<<endl;		
    
    
	if (string_length_%2) {
	    int a_mid = (string_length_-1)/2;
	    
		// odd n_j
		for (int a=0; a < a_mid; ++a) {
			new_delta[a]=0;
		    for (int b=a; b < a_mid; ++b) {
				new_delta[a] += del_term[b];
			}
		}
		
		// string center has no deviations.
		new_delta[a_mid] = 0.0;
		
		// set the conjugates
		for (int a=a_mid+1; a<string_length_;++a) {
			new_delta[a] = - new_delta[string_length_-a-1];
		}
	}
	else {
		// even n_j. 
		int a_mid = string_length_/2 - 1;
	    	
		// first, the regular terms.
		for (int a=0; a<a_mid; ++a) {
			new_delta[a]=0;
		    for (int b=a; b<a_mid; ++b) {
				new_delta[a] += del_term[b];
			}
		}

        // apply the inner-sign term.
        bool inner_narrow = innerPairIsNarrow(sum_jx2_, string_length_, number_roots_odd);
		for (int a=0; a<a_mid+1; ++a) {
			// log_rsq should be sum b=1...n/2 log_rsq[b] now
			new_delta[a] += 0.5L*exp(0.5L*log_rsq) * (inner_narrow?-1.L:1.L); 
		}
		
		// set the conjugates
		for (int a=a_mid+1; a<string_length_; ++a) {
			new_delta[a] = - new_delta[string_length_-a-1];
		}
	}
	return true;
}





vector<complex<double> > CentralString::getCleanerJ (
                            const vector<SymRoots*>& all, const int alpha) const
{
    // outputs
	vector<complex<double> > result(string_length_);
	
	for (int a=0; a<string_length_; ++a) {
		double log_rsq = 0.0;
	    int theta = 0;
	
		double y = 0.5*(string_length_-1) - a + delta_[a];
        double del = delta_[a];
	    int pos = (string_length_-1) - 2*a;
	    
		// string self scattering 		
		for (int b=0; b < string_length_; ++b) {
			
	        double y_diff = (delta_[a]-delta_[b]) - (a-b);
	    
			if (b!=a) {    
			    theta -=  (abs(y_diff)>1.) ? sgn(y_diff) : 0;
			    // split position and delta differences for precision in the case where (1-a+b)==0 etc.
    	        log_rsq += log(sq( ((double)(1-a+b) + (delta_[a]-delta_[b])) / ((double)(1+a-b) - (delta_[a]-delta_[b])) ));
            }   
        }		
		
		// scattering with others
		for (int beta=0; beta<all.size(); ++beta) {
		    // self scattering covered above
		    if (beta==alpha) continue;
		
    		//SymRoots* roots_beta = static_cast<SymRoots*>(all[beta]);
            vector<double> other_x = all[beta]->getRealPairs();
            //vector<double> other_y = all[beta]->getImagPairs();
            vector<double> other_del = all[beta]->getImagPairsDelta();
            vector<int> other_pos = all[beta]->getImagPairsPos();
            vector<complex<double> > other_z = all[beta]->getComplexQuartets();
            
            if (all[beta]->hasOrigin()) {
                theta -=  (abs(y)>1.) ? sgn(y) : 0;
                log_rsq +=  log(sq( (1.+y)/(1.-y) ));
            }
            
            for (int k=0; k<other_x.size(); ++k) {
                double x_k = other_x[k];
                
                // +x
				log_rsq +=  log( (sq(-x_k) + sq(1.+y)) / (sq(-x_k) + sq(1.-y)) );
	            // -x
				log_rsq +=  log( (sq(+x_k) + sq(1.+y)) / (sq(+x_k) + sq(1.-y)) );
	        }

            for (int k=0; k<other_del.size(); ++k) {
                
                // split position and delta for precision in case positions are equal
                long double del_dif = del-other_del[k];
                long double del_sum = del+other_del[k];
                
                long double pos_dif = 0.5L*(pos-other_pos[k]);
                long double pos_sum = 0.5L*(pos+other_pos[k]);
                
                // +iy
			    theta -=  (abs(pos_dif+del_dif)>1.) ? sgn(pos_dif+del_dif) : 0;
				log_rsq +=  4.L*atih(pos_dif+del_dif);
	            
	            // -iy
			    theta -=  (abs(pos_sum+del_sum)>1.) ? sgn(pos_sum+del_sum) : 0;
				log_rsq +=  4.L*atih(pos_sum+del_sum);
            }
            /*                    
            for (int k=0; k<other_y.size(); ++k) {
                double y_k = other_y[k];
                // +iy
			    theta -=  (abs(y-y_k)>1.) ? sgn(y-y_k) : 0;
				log_rsq +=  log( sq(1. + (y-y_k)) / sq(1. - (y-y_k)) );
	            // -iy
			    theta -=  (abs(y+y_k)>1.) ? sgn(y+y_k) : 0;
				log_rsq +=  log( sq(1. + (y+y_k)) / sq(1. - (y+y_k)) );
	        }
            */
                   
            for (int k=0; k<other_z.size(); ++k) {
                double x_k = real(other_z[k]);
                double y_k = imag(other_z[k]);
                // x+iy
				log_rsq +=  log( sq(-x_k) + sq(1. + (y-y_k)) ) - log( sq(-x_k) + sq(1. - (y-y_k)) );
	            // -x-iy
				log_rsq +=  log( sq(+x_k) + sq(1. + (y+y_k)) ) - log( sq(+x_k) + sq(1. - (y+y_k)) );
                // x-iy
				log_rsq +=  log( sq(-x_k) + sq(1. + (y+y_k)) ) - log( sq(-x_k) + sq(1. - (y+y_k)) );
	            // -x+iy
				log_rsq +=  log( sq(+x_k) + sq(1. + (y-y_k)) ) - log( sq(+x_k) + sq(1. - (y-y_k)) );
	        }
	    }
        
        
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
    vector< complex<double> > dirty_j = getCleanerJ(all, alpha);
    vector<int> result(dirty_j.size(), 0);
    sq_error=0.;
    for (int i=0; i<result.size(); ++i) {
        result[i] = round(2.*real(dirty_j[i]));
        sq_error += norm(1.*result[i] - 2.*dirty_j[i])/(double)dirty_j.size();           
    }      
    return result;
}
