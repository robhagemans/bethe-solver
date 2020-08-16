#include "symmetric-string.h"

using namespace std;



class KiteString : public SymRoots {
public:

    KiteString(const int big_n, const int string_length);// const int sum_jx2 );
    // solving policy, default is used by constructor
    void setPolicy(const double initial_deviation=1e-10, const double damping_delta=0.);
    
    virtual std::vector<CVal> getRoots() const;
    virtual int size() const;
    
    virtual std::vector<CVal> getComplexQuartets() const;
    virtual std::vector<double> getRealPairs() const; 

    virtual std::vector<IVal> getImagPairs() const; 
        
    virtual IterResult iterate(const std::vector<SymRoots*>& all, const int alpha, double& convergence);
    virtual IterResult initiate(const std::vector<SymRoots*>& all, const int alpha);
    virtual void refresh();


private:
    double stepRapidity (const std::vector<SymRoots*>& all, const int alpha);

    bool stepDeviation (const std::vector<SymRoots*>& all, const int alpha, 
                        std::vector<double>& new_delta, std::vector<double>& new_epsilon);

    void setInitialValues(const double initial_eps, const double initial_del);

    // chain length
    double big_n_;
    
    // number of roots in this string
    int string_length_;          
    
    // sum of bethe quantum numbers for the string, times 2
    //int sum_jx2_; 

    bool number_roots_odd_;
    
    std::vector<double> epsilon_;   
    std::vector<double> delta_;   

    std::vector<double> new_epsilon_;   
    std::vector<double> new_delta_;   
    
    std::vector<double> last_epsilon_;   
    std::vector<double> last_delta_;   
    
    std::vector<double> step_epsilon_;   
    std::vector<double> step_delta_;   
    
    
    // hold a deviation - do not update
    std::vector<bool> hold_;


    // number of iteration steps passed    
    int iterations_;
    
    /* solving policy */

    // size of small deviation to initialize
    double initial_deviation_;
    // maximum number of steps spent in 'running' phase (exponential growing of lambda to find solvable deviation) 
    //int run_max_; 
    // damping factor for deviance delta
	double damping_delta_; 
};




/** class KiteString **/


KiteString::KiteString(const int big_n, const int string_length)//, const int sum_jx2)
      : big_n_(big_n), //sum_jx2_(sum_jx2), 
        string_length_(string_length),  number_roots_odd_(false),
        epsilon_(string_length), delta_(string_length), 
        new_epsilon_(string_length), new_delta_(string_length), 
        last_epsilon_(string_length), last_delta_(string_length), 
        step_epsilon_(string_length), step_delta_(string_length), 
        hold_(string_length, false),
        iterations_(0) 
{
    setPolicy();
}

void KiteString::setPolicy(const double initial_deviation, const double damping_delta)
{
    initial_deviation_ = initial_deviation;
    damping_delta_ = damping_delta;
}
    

vector<CVal> KiteString::getComplexQuartets() const
{
    vector<CVal> quartets(string_length_/2-1);
    // run through half string, excluding odd strings' centres which are real.
    // excluding top & bootom which are imag.
    for (int a=1; a<string_length_/2; ++a) {
        //quartets[a-1] = epsilon_[a] + I*(delta_[a] + 0.5*(string_length_-1) - a);
        quartets[a-1] = { 0., epsilon_[a], (string_length_-1) - 2*a,  delta_[a] };
    }
    return quartets;
}

vector< double > KiteString::getRealPairs() const
{
    // no real root for even string
    if (string_length_%2) return { epsilon_[string_length_/2] };
    else return { };
}


vector<IVal> KiteString::getImagPairs() const
//{   return { delta_[0] + 0.5*(string_length_-1) };  }
{   return { {(string_length_-1), delta_[0]} };  }

    
vector<CVal> KiteString::getRoots() const
{
    vector<CVal> roots (string_length_*2-2);
    roots[0] = { 0., 0., (string_length_-1), delta_[0] };
    roots[string_length_*2-3] =  { 0., 0., -(string_length_-1), -delta_[0] };
    
    for (int a=1; a<string_length_-1; ++a) {
        // run through string, this covers conjugates.
        roots[a] = 
                { 0., epsilon_[a], (string_length_-1) - 2*a, delta_[a]};
        // opposite string
        roots[string_length_*2-3-a] = 
                { 0., -epsilon_[a], -(string_length_-1) + 2*a, -delta_[a]};
    }
    return roots;
}


int KiteString::size() const
{   return string_length_*2-2;  }




/** calculate a new value for the deviance (delta) and aberration (epsilon) **/
bool KiteString::stepDeviation (
                const vector<SymRoots*>& all, const int alpha, 
                vector<double>& new_delta, vector<double>& new_epsilon)
{
	if (string_length_<3) return true;

	long double theta_sum = 0.0;
	long double log_rsq_sum = 0.0;

	const int a_mid = (string_length_-1)/2;
	vector<long double> del_term (string_length_, 0.);
	vector<long double> eps_term (string_length_, 0.);
	
	vector<long double> theta_store (string_length_, 0.);
	vector<long double> log_rsq_store (string_length_, 0.);
	
	for (int a=0; a<string_length_; ++a) {
		 
	    long double theta = 0.0;
	    long double log_rsq = 0.0;
		 
		const CVal z = { 0., epsilon_[a], (string_length_-1) - 2*a, delta_[a] };  
		const double x = epsilon_[a];
	    //const double y = 0.5*(string_length_-1) - a + delta_[a];
	    //const double pos = 0.5*(string_length_-1) - a;
	    const int posx2 = (string_length_-1) - 2*a;
	    const double del = delta_[a];
	    
		// string self scattering minus dangerous terms
		// exclude last root as it's covered by the opposite of the first		
		for (int b=0; b < string_length_-1; ++b) {
		    if (b==a) continue;
		    // take care to separate position and delta as we're keeping 'dangerous' scattering with nearby opposite string
			
		    const double x_sum = epsilon_[a] + epsilon_[b];
            const int pos_sum = (string_length_-1) - (a+b);
		    const double del_sum = delta_[a] + delta_[b];

            const double x_diff = epsilon_[a] - epsilon_[b];
            const int pos_diff = -(a-b);
            const double del_diff = (delta_[a]-delta_[b]);

		    // take care to separate pos and delta as we're keeping 'dangerous' scattering with nearby opposite string
		    if (a-b!=1)  theta -=  xi( x_diff, (1+pos_diff) + del_diff );
		    if (a-b!=-1) theta -=  xi( x_diff, (1-pos_diff) - del_diff );
		    if (a-b!=1)  log_rsq +=  log( sq(x_diff) + sq((1+pos_diff) + del_diff) );
		    if (a-b!=-1) log_rsq += -log( sq(x_diff) + sq((1-pos_diff) - del_diff) );
            
	        // opposite string, excluding opposite root absorbed in kinetic term
		    theta -=  xi( x_sum, (1+pos_sum)+del_sum )  +  xi( x_sum, (1-pos_sum)-del_sum );
		    log_rsq +=  log(  (sq(x_sum) + sq((1+pos_sum)+del_sum))  /  (sq(x_sum) + sq((1-pos_sum)-del_sum))  );
		}		

	    for (int beta=0; beta<all.size(); ++beta) {
	        // exclude self scattering
	        if (beta==alpha) continue;
	    	// scattering with others
	    	bool other_o = all[beta]->hasOrigin();
	    	vector<double> other_x = all[beta]->getRealPairs();
	    	vector<IVal> other_y = all[beta]->getImagPairs();
	    	vector<CVal> other_z = all[beta]->getComplexQuartets();
            
            theta -= sum_re_atan_x2(z, other_o, other_x, other_y, other_z);
            log_rsq += sum_im_atan_x4(z, other_o, other_x, other_y, other_z);
        }
                
		// contribution M mod 2 from sum of conjugate quantum numbers
		theta -= PI*number_roots_odd_; 
		
		// kinetic theta, absorbing opposite root
		theta += (big_n_-1.) * (xi(2.*x, (1+posx2)+2.*del)  +  xi(2.*x, (1-posx2)-2.*del));
		log_rsq -= (big_n_-1.) * log( (sq(2.*x) + sq((1+posx2)+2.*del))  /  (sq(2.*x) + sq((1-posx2)-2.*del)) );


        theta_sum += theta;
        log_rsq_sum += log_rsq;
		// theta is now sum b=1..a theta[b]
		// log_r is now sum b=1..a log_r[b]
		del_term[a] = -cos(theta_sum)*exp(0.5*log_rsq_sum);
		eps_term[a] = sin(theta_sum)*exp(0.5*log_rsq_sum);
        
        theta_store[a] = theta;
        log_rsq_store[a] = log_rsq;
	}	

	for (int a=string_length_-2; a>=0; --a) {   
	    theta_store[a] += theta_store[a+1];
	    log_rsq_store[a] += log_rsq_store[a+1];
	}
		
    // outputs
    new_delta.assign(string_length_, 0.);
    new_epsilon.assign(string_length_, 0.);

    // add differences
    //new_epsilon[0] = 0.;
    for (int a=1; a<=a_mid; ++a) {   
	    new_epsilon[a] += new_epsilon[a-1] - eps_term[a-1];
    }
    //new_delta[a_mid] = 0.;
    for (int a=a_mid-1; a>=0; --a) {   
	    new_delta[a] += new_delta[a+1] + del_term[a];
	}

/*
	if (0==string_length_%2) {
		// even n_j: apply the inner-sign term.
        bool inner_narrow = innerPairIsNarrow(sum_jx2_, string_length_, number_roots_odd_);
		for (int a=0; a<a_mid+1; ++a) {
			// log_r should be sum b=1...n/2 log_r[b] now
			new_delta[a] += 0.5*exp(0.5*log_rsq) * (inner_narrow?-1.:1.); 
		}
	}
*/

	// set the conjugates
	for (int a=a_mid+1; a<string_length_-1; ++a) {
		new_delta[a] = - new_delta[string_length_-a-1];
		new_epsilon[a] = new_epsilon[string_length_-a-1];
	}


	return true;
}





IterResult KiteString::iterate(const vector<SymRoots*>& all, const int alpha, double& convergence)
{
	vector<double> step_delta;
	vector<double> step_epsilon;

	
    const int a_mid = (string_length_-1)/2; 

	double eps0 = stepRapidity(all, alpha);
 
    if (!stepDeviation(all, alpha, step_delta, step_epsilon)) 
        return IterResult::iter_err;

    
    for (int a=1; a<=a_mid;++a) {
        step_epsilon[a] += eps0 - step_epsilon[a_mid];
    }
    
    if (!string_length_%2) step_epsilon[a_mid+1] = step_epsilon[a_mid];
    step_delta[string_length_-1] = -step_delta[0];

    
    for (int a=0; a<a_mid;++a) {
        new_delta_[a] = dampedStep(delta_[a], step_delta[a], damping_delta_);
                        //secantStep(last_delta_[a], delta_[a], step_delta_[a], step_delta[a]);
        
    }

    for (int a=1; a<=a_mid;++a) {
        new_epsilon_[a] = dampedStep(epsilon_[a], step_epsilon[a], damping_delta_);
                          //secantStep(last_epsilon_[a], epsilon_[a], step_epsilon_[a], step_epsilon[a]);
    }
    new_delta_[a_mid] = 0.;
    for (int a=0; a<=a_mid; ++a) {
        new_delta_[string_length_-1-a] = -new_delta_[a];
        new_epsilon_[string_length_-1-a] = new_epsilon_[a];
    }   
    
    step_delta_ = step_delta;
    step_epsilon_=step_epsilon; 

    // enforce holds
    for (int a=0; a<string_length_;++a) 
        if (hold_[a]) {
		    new_delta_[a] = delta_[a];
		    new_epsilon_[a] = epsilon_[a];
	    }
    
    // calculate convergence
	for (int a=0; a < string_length_; ++a)
		convergence += sq(new_delta_[a] - delta_[a]) + sq(new_epsilon_[a] - epsilon_[a]);

    ++iterations_;
	return IterResult::iter_ok;
}



void KiteString::refresh()
{
	delta_ = new_delta_;
	epsilon_ = new_epsilon_;
    return;
}

void KiteString::setInitialValues(const double initial_eps, const double initial_del)
{
    const int a_mid = (string_length_-1)/2;
    for(int a=0; a<a_mid; ++a) {
        delta_[a] = initial_del * (string_length_ - 2*a);
        delta_[string_length_-1-a] = -delta_[a];
    }
    for(int a=1; a<=a_mid; ++a) {
        epsilon_[a] = epsilon_[string_length_-1-a] = initial_eps;
    }
    epsilon_[0] = epsilon_[string_length_-1] = 0;    	
    	
    // for odd strings, set central deviations.
    // for even strings, ensure the deviation of the inner pair has the correct sign
    if (string_length_%2) {
        delta_[a_mid] = 0.;
    }
    else {
/*        if (innerPairIsNarrow(sum_jx2_, string_length_, number_roots_odd_)) {
            delta_[a_mid] = -initial_deviation_;
            delta_[a_mid+1] = -delta_[a_mid];
        } 
*/
    }

}

	
IterResult KiteString::initiate(const vector<SymRoots*>& all, const int alpha)
{
    number_roots_odd_ = false;
    for (int beta=0; beta<all.size();++beta) { 
        if (all[beta]->hasOrigin()) {
            number_roots_odd_=true;
            break;
        } 
    }
    //setInitialValues(0.01, 0.501);
    setInitialValues(initial_deviation_, initial_deviation_/10.);
    if (!stepDeviation(all, alpha, step_delta_, step_epsilon_)) 
        return IterResult::iter_err;

    if (!string_length_%2) step_epsilon_[(string_length_+1)/2] = step_epsilon_[(string_length_-1)/2];
    step_delta_[string_length_-1] = -step_delta_[0];
    
    last_delta_ = delta_;
    last_epsilon_ = epsilon_;
    //setInitialValues(0.02, 0.502);
    setInitialValues(2.*initial_deviation_, 2.*initial_deviation_/10.);
    
    return IterResult::iter_ok;
}    


/** calculate a new value for the string centre (lambda) **/
double KiteString::stepRapidity (const vector<SymRoots*>& all, const int alpha)
{
    const double lambda_ = 0;
    const int sum_jx2_ = -1;
    
    
	long double scatter = 0.0;		
	for (int a=0; a < string_length_; ++a) {

		const CVal z = { 0., epsilon_[a], (string_length_-1) - 2*a, delta_[a] };  
		
        const double x = lambda_ + epsilon_[a];
        const double y = 0.5*(string_length_-1) - a + delta_[a];
            
        // scattering with opposite strings
	    for (int b=1; b < string_length_-1; ++b) {
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
	    	const vector<IVal> other_y = all[beta]->getImagPairs();
	    	const vector<CVal> other_z = all[beta]->getComplexQuartets();
            scatter += sum_re_atan_x2( z, other_o, other_x, other_y, other_z );
        }
        
	}
	// 2\pi \sum_{a=1}^n J_a + \theta\sub{other}^a = \sum_{a=1}^n \theta\sub{kin}^a 
    // absorbed opposite root in N-1 for faster convergence
    const long double phase = 0.5*(PI*sum_jx2_ + scatter) / (big_n_-1.);
    const int a_mid = (string_length_-1)/2;

	long double sum_kin = 0.0;
    for (int a=0; a< a_mid; ++a) {
        sum_kin +=  xi(2.*(lambda_ + epsilon_[a]), 1.+ (0.5*(string_length_-1) - a + delta_[a]) )
                    + xi(2.*(lambda_ + epsilon_[a]), 1.-( 0.5*(string_length_-1) - a + delta_[a]) );

	}
	// missing term: atan (2*x);
    return 0.5*tan( phase-sum_kin ) ;
}









