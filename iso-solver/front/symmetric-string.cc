#include "symmetric-string.h"

using namespace std;




// determine width of inner pair
// only works for even N
const bool innerPairIsNarrow(const int sum_jx2, const int string_length, const int number_roots_odd)
{
    // the question doesn't make sense for odd strings.
    if (string_length%2) return false;

    // inner pair is wide if the bethe quantum numbers are J, J+1
    // if M odd, all Js are integer, so 2J is even, so 2J+1 is odd. 
    // if M even, all Js are half-odd-integer, so 2J is odd, so 2J+1 is even.
    // so the parities are the same:  2J+1 == M (mod 2)
    // every other pair in the string is wide, a 2J_n+1 is added for each
    // for an l-string with l even, there are l/2 pairs so l/2-1 other pairs.
    // each other pair adds M (mod 2).
    // 0.5*sum_2J == 2J+1 - [string is narrow] + sum_{other n} 2J_n+1 
    // 0.5*sum_2J == M - [string is narrow] + (l/2-1)*M   (mod 2)
    // [string is narrow] =  M*(l/2) - sum_2J         (mod 2)  
    return ( number_roots_odd*(string_length/2) - sum_jx2/2 )%2 ;    
}





/** class String **/


SymString::SymString(const int big_n, const int string_length, const int sum_jx2, const double lambda)
      : big_n_(big_n), 
        string_length_(string_length), sum_jx2_(sum_jx2),  number_roots_odd_(false),
        lambda_(lambda), epsilon_(string_length), delta_(string_length), 
        new_lambda_(0.), new_epsilon_(string_length), new_delta_(string_length), 
        hold_(string_length, false), kite_mate_(0), mate_beta_(-1), n_on_axis_(0),
        iterations_(0) 
{
    setPolicy();
}

void SymString::setPolicy(const double initial_deviation, const int steps_no_deviation, const double damping_delta)
{
    initial_deviation_ = initial_deviation;
    steps_no_deviation_ = steps_no_deviation;
    damping_delta_ = damping_delta;
}
    
bool SymString::setOnAxis(const int n_on_axis)
{
    n_on_axis_ = n_on_axis;
    return n_on_axis < (string_length_+1)/2;
}    


bool SymString::coupleKite(Kite& mate, const int mate_beta)
{
    // set mate reference
    kite_mate_ = &mate;
    mate_beta_ = mate_beta;
    return string_length_%2; // even strings not impl yet
} 


vector<CVal> SymString::getComplexQuartets() const
{
    vector<CVal> result(string_length_/2);
    // run through half string, excluding odd strings' centres which are real.
    for (int a=n_on_axis_; a<string_length_/2; ++a) {
        //result[a] = lambda_ + epsilon_[a] + I*(delta_[a] + 0.5*(string_length_-1) - a);
        result[a] = { lambda_, epsilon_[a], (string_length_-1) - 2*a, delta_[a] };
    
    }
    return result;
}

vector< double > SymString::getRealPairs() const
{
    // no real root for even string
    if (!kite_mate_&& string_length_%2) return { lambda_ };
    else return { };
}


vector<IVal> SymString::getImagPairs() const
{
    vector<IVal> result(n_on_axis_);
    // run through half string, excluding odd strings' centres which are real.
    for (int a=0; a<n_on_axis_; ++a) {
        result[a] = { (string_length_-1) - 2*a, delta_[a] };
    }
    return result;
}


vector<CVal> SymString::getRoots() const
{
    const int vec_length = (string_length_-n_on_axis_)*2 - ((kite_mate_)? 2 : 0);
    const int left_length = string_length_ - ((kite_mate_)? 1 : 0);
    vector<CVal> roots (vec_length);
    for (int a=0; a<n_on_axis_; ++a) {
        // top  // 0
        roots[a] =                          { 0.,0.,  (string_length_-1) - 2*a,  delta_[a] };
        // bottom // 3
        roots[left_length-1-a] =            { 0.,0., -(string_length_-1) + 2*a, -delta_[a] };
    }
    const int a_last = (kite_mate_) ? (string_length_)/2 : (string_length_+1)/2;
    for (int a=n_on_axis_; a<a_last; ++a) {
        // top left // 1
        roots[a] =                          { lambda_,  epsilon_[a],  (string_length_-1) - 2*a,  delta_[a] };
        // bottom left // 2
        roots[left_length-1-a] =            { lambda_,  epsilon_[a], -(string_length_-1) + 2*a, -delta_[a] };
        // top left // 4
        roots[left_length-n_on_axis_+a] =   {-lambda_, -epsilon_[a],  (string_length_-1) - 2*a,  delta_[a] };
        // bottom left // 5
        roots[vec_length-1-a] =             {-lambda_, -epsilon_[a], -(string_length_-1) + 2*a, -delta_[a] };
    }
    return roots;
}


int SymString::size() const
{
    return 2*(string_length_ - n_on_axis_ - (kite_mate_?1:0));
}



/** calculate a new value for the deviance (delta) and aberration (epsilon) **/
bool SymString::stepDeviation (
                const vector<SymRoots*>& all, const int alpha, 
                vector<double>& new_delta, vector<double>& new_epsilon)
{
	// real rapidity, no deviations.
	if (1==string_length_) 		{
	    new_delta[0] = new_epsilon[0] = { 0. };
		return true;
	}

	HPIVal theta = { 0, 0. };
	double log_rsq = 0.0;
	vector<double> del_term (string_length_/2, 0.);
	vector<double> eps_term (string_length_/2, 0.);

    // even: runs to a_mid+1. odd: runs to a_mid.
	for (int a=0; a<string_length_/2; ++a) {
		
		const CVal z = {lambda_, epsilon_[a], (string_length_-1)-2*a, delta_[a]};
		
		for (int b=0; b < string_length_-n_on_axis_; ++b) {
			if (b==a) continue;
			// off-axis (b==n_on_axis..l-n_on_axis) : both upper and lower root of each pair included in loop
			// on-axis (b==0..n_on_axis-1) : only upper root included in loop, lower root will be opposite
			const CVal z_b = {lambda_, epsilon_[b], (string_length_-1)-2*b, delta_[b]};
            
            if (a-b!=1)  {
                theta.sub(xi_plus( dif(z, z_b) ));
		        log_rsq +=  zeta_plus_x2( dif(z, z_b) );
		    }
		    if (a-b!=-1) {
		        theta.sub(xi_minus( dif(z, z_b) ));
		        log_rsq += -zeta_minus_x2( dif(z, z_b) );
            }
	        // opposite string, excluding opposite root absorbed in kinetic term
            theta.sub(mul( re_atan( sum(z, z_b) ), 2));
            log_rsq += im_atan_x4( sum(z, z_b) );
		}		

	    for (int beta=0; beta<all.size(); ++beta) {
	        // exclude self scattering
	        if (beta==alpha) continue;
            // exclude scattering with mate
	    	if (kite_mate_ && beta==mate_beta_) continue;
	    	// scattering with others
	    	const bool other_o = all[beta]->hasOrigin();
	    	const vector<double> other_x = all[beta]->getRealPairs();
	    	const vector<IVal> other_y = all[beta]->getImagPairs();
	    	const vector<CVal> other_z = all[beta]->getComplexQuartets();
	    	theta.sub(mul( sum_re_atan( z, other_o, other_x, other_y, other_z ), 2));
            log_rsq += sum_im_atan_x4( z, other_o, other_x, other_y, other_z );
        }
        

        if (kite_mate_) {
            // only include top and bottom here, we scatter with the centre through ourselves.
             
            // top kite root, excluding dangerous terms
	        const IVal y_top = { kite_mate_->pos_, kite_mate_->del_ };
	        theta.sub( mul(re_atan( dif(z, y_top) ),2));
	        theta.sub( mul(re_atan( sum(z, y_top) ),2));
	        log_rsq += im_atan_x4( dif(z, y_top) ) + im_atan_x4( sum(z, y_top) );
        }
        
		// contribution M mod 2 from sum of conjugate quantum numbers
		//theta -= PI*number_roots_odd_; 
		theta.sub_half_pi( 2*number_roots_odd_ ); 

		// kinetic theta, absorbing opposite root
		//theta += (big_n_-1)* re_atan_x2( mul(z, 2) );
		theta.add(mul( re_atan( mul(z, 2) ), 2*(big_n_-1)));
		log_rsq -= (big_n_-1)* im_atan_x4( mul(z, 2) );

		// theta is now sum b=1..a theta[b]
		// log_r is now sum b=1..a log_r[b]
		del_term[a] = -cos(theta)*exp(0.5*log_rsq);
		eps_term[a] = sin(theta)*exp(0.5*log_rsq);
	}	
		
    // outputs
	const int a_mid = (string_length_-1)/2;
    new_delta.assign(string_length_, 0.);
    new_epsilon.assign(string_length_, 0.);

    // add differences
    for (int a=a_mid-1; a>=0; --a) {   
	    new_delta[a] += new_delta[a+1] + del_term[a];
        //if (a>=n_on_axis_) 
        new_epsilon[a] += new_epsilon[a+1] + eps_term[a];
    }

	//new_epsilon[a_mid] = 0.0;
	//new_delta[a_mid] = 0.0;
	if (0==string_length_%2) {
		// even n_j: apply the inner-sign term.
        const bool inner_narrow = innerPairIsNarrow(sum_jx2_, string_length_, number_roots_odd_);
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
double SymString::stepRapidity (const vector<SymRoots*>& all, const int alpha)
{
    // kite determines rapidity if set 
    if (kite_mate_) return real(kite_mate_->z_next_);

	HPIVal scatter { 0, 0. };		
	for (int a=0; a < string_length_; ++a) {

		const CVal z = {lambda_, epsilon_[a], (string_length_-1)-2*a, delta_[a]};

        // scattering with opposite strings
	    for (int b=0; b < string_length_; ++b) {
            // direct opposite - absorb into kinetic term
            if (b+a==string_length_-1) continue;
            const CVal z_b = {-lambda_, -epsilon_[b], (string_length_-1) - 2*b, delta_[b] };
 			scatter.add(mul(re_atan( dif(z, z_b) ), 2));
		}

        // scattering with others
        for (int beta=0; beta<all.size(); ++beta) {
	        // exclude self scattering
	        if (beta==alpha) continue;
	    	const bool other_o = all[beta]->hasOrigin();
	    	const vector<double> other_x = all[beta]->getRealPairs();
	    	const vector<IVal> other_y = all[beta]->getImagPairs();
	    	const vector<CVal> other_z = all[beta]->getComplexQuartets();
            //scatter.add(sum_re_atan_x2( z, other_o, other_x, other_y, other_z ));
            scatter.add(mul(sum_re_atan( z, other_o, other_x, other_y, other_z ), 2));
        }
	}
	
/*
    const int a_mid = (string_length_-1)/2;
    double sum_kin = 0.0;
	
	// term for even strings
	if (0==string_length_%2)   sum_kin += xi_CAUTION( lambda_, -delta_[a_mid] ); 
	
    for (int a=0; a< a_mid; ++a) {
		sum_kin += 	  xi_CAUTION( lambda_ + epsilon_[a+1], 0.5*string_length_ - (a+1) + delta_[a+1]  )
					+ xi_CAUTION( lambda_ + epsilon_[a],  -0.5*string_length_ + (a+1) - delta_[a]  );
	}


	// 2\pi \sum_{a=1}^n J_a + \theta\sub{other}^a = \sum_{a=1}^n \theta\sub{kin}^a 
    // absorbed opposite root in N-1 for faster convergence
    const double phase_x2 = -sum_kin_x2.val() + (PI*sum_jx2_ + scatter.val()) / (big_n_-1.);
    return -epsilon_[0] + ( 0.5*string_length_ + delta_[0] )*tan( 0.5*phase_x2 ) ;
*/

    HPIVal sum_kin_x2 = { 0, 0. };
    for (int a=0; a < string_length_; ++a) {
        const CVal z = { lambda_, epsilon_[a], (string_length_-1) - 2*a, delta_[a] };
	    if (string_length_-1!=a) sum_kin_x2.add( xi_minus( mul(z,2) ) );
	    if (0!=a)  sum_kin_x2.add( xi_plus( mul(z,2) ) );
    }
    
	// 2\pi \sum_{a=1}^n J_a + \theta\sub{other}^a = \sum_{a=1}^n \theta\sub{kin}^a 
    // absorbed opposite root in N-1 for faster convergence
    scatter.add_half_pi(2*sum_jx2_);
    scatter.sub(mul( sum_kin_x2, big_n_-1 ));
    
    return -epsilon_[0] + ( 0.5*string_length_ + delta_[0] )*tan(div(scatter, 2*(big_n_-1))) ;
}




IterResult SymString::iterate(const vector<SymRoots*>& all, const int alpha, double& convergence)
{

//cerr<<lambda_<<" "<<epsilon_<<delta_<<endl;

	new_lambda_ = stepRapidity(all, alpha);
	if (!finite(new_lambda_))  
	    return IterResult::iter_err_nonfinite;

	
    if (iterations_ >= steps_no_deviation_) {
	    if (!stepDeviation(all, alpha, new_delta_, new_epsilon_)) 
	        return IterResult::iter_err;
        for (int a=0; a<string_length_; ++a) 
	        new_delta_[a] = dampedStep(delta_[a], new_delta_[a], damping_delta_);
    }

//cerr<<new_lambda_<<" "<<new_epsilon_<<new_delta_<<endl;
    if (n_on_axis_) {
        // reconcile lambda and epsilon
        new_lambda_ = 0.5*(new_lambda_-new_epsilon_[n_on_axis_-1]);   
        
        // fix on-axis roots to axis
        for(int a=0; a<n_on_axis_; ++a) {
            new_epsilon_[a] = -new_lambda_;
            new_epsilon_[string_length_-1-a] = -new_lambda_;
        }
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
	return IterResult::iter_ok;
}



void SymString::refresh()
{
	lambda_ = new_lambda_;
	delta_ = new_delta_;
	epsilon_ = new_epsilon_;
    return;
}

	
IterResult SymString::initiate(const vector<SymRoots*>& all, const int alpha)
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
            lambda_ = abs(0.5*sum_jx2_/double((string_length_-2*n_on_axis_)*big_n_));    
    }
    
    const int a_mid = (string_length_-1)/2;
    for(int a=0; a<a_mid; ++a) {
        epsilon_[a] = epsilon_[string_length_-1-a] = initial_deviation_;
        delta_[a] = initial_deviation_;
        delta_[string_length_-1-a] = -delta_[a];
    }
    epsilon_[a_mid] = 0.;

    for(int a=0; a<n_on_axis_; ++a) {
        epsilon_[a] = -lambda_;
        epsilon_[string_length_-1-a] = -lambda_;
    } 
    	
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
    
    return IterResult::iter_ok;
}


