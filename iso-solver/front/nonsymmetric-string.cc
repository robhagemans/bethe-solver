#include "nonsymmetric-string.h"

using namespace std;



/** class String **/


NonsymString::NonsymString(const int big_n, const int string_length, const int ix2)
      : big_n_(big_n),  string_length_(string_length), ix2_(ix2), have_inner_sign_(false),
        number_roots_(0), inner_narrow_(false),
        lambda_(0.), epsilon_(string_length), delta_(string_length),
        new_lambda_(0.), new_epsilon_(string_length), new_delta_(string_length),
        last_lambda_(0.), last_epsilon_(string_length), last_delta_(string_length),
        step_lambda_(0.), step_epsilon_(string_length), step_delta_(string_length),
        hold_(string_length, false),
        iterations_(0),
        initial_deviation_(1e-10),
        damping_dev_(0)
{ }

vector<CVal> NonsymString::getComplexPairs() const
{
    vector<CVal> pairs(string_length_/2);
    // run through half string, excluding odd strings' centres which are real.
    for (int a=0; a<string_length_/2; ++a) {
        pairs[a] = { lambda_, epsilon_[a], (string_length_-1) - 2*a, delta_[a] };
    }
    return pairs;
}


vector< double > NonsymString::getRealRoots() const
{
    // no real root for even string
    if (string_length_%2) return { lambda_ };
    else return {};
}


vector<double> NonsymString::getRapidities() const
{   return { lambda_ };  }

int NonsymString::stringLength() const
{   return string_length_; }

int NonsymString::size() const
{   return string_length_; }


vector<CVal> NonsymString::getRoots() const
{
    vector<CVal> roots (string_length_);
    for (int a=0; a<string_length_; ++a) {
        // run through string, this covers conjugates.
        roots[a] = {lambda_, epsilon_[a], (string_length_-1) - 2*a ,delta_[a] };
    }
    return roots;
}


/** calculate a new value for the deviance (delta) and aberration (epsilon) **/
bool NonsymString::stepDeviation (
                const vector<NonsymRoots*>& all, const int alpha,
                vector<double>& new_delta, vector<double>& new_epsilon)
{
    // real rapidity, no deviations.
	if (1==string_length_) {
		new_delta = new_epsilon = { 0. };
		return true;
	}

	HPIVal theta = { 0, 0. };
	double log_rsq = 0.0;
	vector<double> del_term (string_length_/2, 0.);
	vector<double> eps_term (string_length_/2, 0.);

    // 4-> 2 (0 1 2 3),   3 -> 1 (0 1 2)
	for (int a=0; a<string_length_/2; ++a) {

		const CVal z = { lambda_, epsilon_[a], (string_length_-1) - 2*a, delta_[a] };
		// string self scattering minus dangerous terms
		for (int b=0; b < string_length_; ++b) {
            if (b==a) continue;

		    const CVal z_b = { lambda_, epsilon_[b], (string_length_-1) - 2*b, delta_[b] };

		    if (a-b!=1) {
		        theta.sub(  xi_plus( dif(z, z_b) ));
		        log_rsq +=  zeta_plus_x2( dif(z, z_b) );
            }
		    if (a-b!=-1) {
		        theta.sub(  xi_minus( dif(z, z_b) ));
		        log_rsq += -zeta_minus_x2( dif(z, z_b) );
		    }
		}

		// scattering with others
		for (int beta=0; beta<all.size(); ++beta) {
		    // self scattering covered above
		    if (beta==alpha) continue;
	        const vector<double> other_x = all[beta]->getRealRoots();
            const vector<CVal> other_z = all[beta]->getComplexPairs();
            //theta.sub( sum_re_atan_x2( z, other_x, other_z ) );
            theta.sub(mul( sum_re_atan( z, other_x, other_z ), 2));

            log_rsq += sum_im_atan_x4( z, other_x, other_z );
	    }

		// contribution M mod 2 from sum of conjugate quantum numbers
		//theta -= PI*number_roots_; //%2
		theta.sub_half_pi( 2*number_roots_ );

		// kinetic theta
        theta.add(mul( re_atan( mul(z,2) ), 2*big_n_ ));
		log_rsq -= big_n_ * im_atan_x4( mul(z,2) );

		// theta is now sum b=1..a theta[b]
		// log_rsq is now sum b=1..a log_r[b]
		del_term[a] = -cos(theta)*exp(0.5*log_rsq);
		eps_term[a] = sin(theta)*exp(0.5*log_rsq);
	}

	// outputs
	new_delta.assign(string_length_, 0.);
	new_epsilon.assign(string_length_, 0.);

    // sum the delta and epsilon terms
    // 4-> 1 (0 1 2 3),   3 -> 1 (0 1 2)
    const int a_mid = (string_length_-1)/2;

    // add differences
    for (int a=a_mid-1; a>=0; --a) {
	    new_delta[a] += new_delta[a+1] + del_term[a];
        new_epsilon[a] += new_epsilon[a+1] + eps_term[a];
    }

	//new_epsilon[a_mid] = 0.0;
	//new_delta[a_mid] = 0.0;
	if (0==string_length_%2) {

if (inner_narrow_!=	 innerPairIsNarrow(all, alpha))
cerr<<"inner sign change"<<endl;

	    // recalculate inner sign
        inner_narrow_ = innerPairIsNarrow(all, alpha);
		// even n_j: apply the inner-sign term.
		for (int a=0; a<a_mid+1; ++a) {
			// log_r should be sum b=1...n/2 log_r[b] now
			new_delta[a] += 0.5*exp(0.5*log_rsq) * (inner_narrow_?-1.:1.);
		}
	}

	// set the conjugates
	for (int a=a_mid+1; a<string_length_; ++a) {
		new_delta[a] = - new_delta[string_length_-a-1];
		new_epsilon[a] = new_epsilon[string_length_-a-1];
	}

	// signal success
	return true;
}





inline HPIVal re_atan_BT_plus (const CVal& z)
{
    if (z.pos==2) {
        if (z.del==0.) return { 0, 0 };
//          else return re_atan(z.real()/(0.5*(2-z.pos)-z.del)); // close to zero if lam+eps is close to 0, ie assumption eps<<del.

        else return re_atan(-(0.5*(2-z.pos)-z.del)/z.real()); // close to zero if delta is close to 0, ie assumption del<<eps.
    }
    return re_atan(z.real()/(0.5*(2-z.pos)-z.del));
}

HPIVal sum_re_atan_BT_plus (const CVal& z, const vector<double>& other_x, const vector<CVal>& other_z)
{
    HPIVal theta = { 0, 0. };
    for (int k=0; k<other_x.size(); ++k) {
        theta.add(re_atan_BT_plus( dif(z, other_x[k]) ));
    }
    for (int k=0; k<other_z.size(); ++k) {
        theta.add(re_atan_BT_plus( dif(z, other_z[k]) ));
		theta.add(re_atan_BT_plus( dif(z, conj(other_z[k])) ));
    }
    return theta;
}

// for zero deviations, this becomes the bethe-takahashi equation
double NonsymString::stepRapidity (const vector<NonsymRoots*>& all, const int alpha)
{
	HPIVal scatter = { 0, 0. };
	for (int a=0; a < string_length_; ++a) {
        const CVal z = { lambda_, epsilon_[a], (string_length_-1) - 2*a, delta_[a] };
		// scattering with others
		for (int beta=0; beta<all.size(); ++beta) {
		    if (beta==alpha) continue;

	        const vector<double> other_x = all[beta]->getRealRoots();
            const vector<CVal> other_z = all[beta]->getComplexPairs();
		    scatter.add(mul( sum_re_atan_BT_plus( z, other_x, other_z ), 2));
		}
	}

    // kinetic term - terms that cancel in B-T approximation only
    // excluded here are the ln-terms with highest imaginary part i(n_j + delta)
    // (that's the string_length-1 term)
    HPIVal sum_kin_x2 = { 0, 0. };
    for (int a=0; a < string_length_-1; ++a) {
        const CVal z = { lambda_, epsilon_[a], (string_length_-1) - 2*a, delta_[a] };

        sum_kin_x2.add( mul(re_atan_BT_plus(mul(z,2)), 2) );
	}

    // 2\pi \sum_{a=1}^n J_a + \theta\sub{other}^a = \sum_{a=1}^n \theta\sub{kin}^a
    scatter.add_half_pi(2*ix2_);
    scatter.sub(mul( sum_kin_x2, big_n_));
    return -epsilon_[0] + ( 0.5*string_length_ + delta_[0] )*tan(div(scatter, 2*big_n_)) ;
}



// determine width of inner pair
// only works for even N
bool NonsymString::innerPairIsNarrow( const std::vector<NonsymRoots*>& all, const int alpha) const
{
    // note that the question doesn't make sense for odd strings, answer undefined if string_length odd.
    int number_odd_shorter = 0;
	int sign_term = 0;
	for (int k=0; k < all.size(); ++k) {

		const int l_k = all[k]->stringLength();
        if (l_k%2 && l_k<string_length_) number_odd_shorter += all[k]->size()/l_k;

		if (k!=alpha && l_k==string_length_) {
            const vector<double> lambda_k = all[k]->getRapidities();
            for (int beta=0; beta<lambda_k.size(); ++beta) {
                sign_term += isgn(lambda_ - lambda_k[beta]);
	        }
	    }
	}

    return ( (ix2_-big_n_-sign_term + string_length_*(big_n_-2*number_roots_-2+string_length_) )/2 + number_odd_shorter )%2 ;
}




IterResult NonsymString::iterate(const vector<NonsymRoots*>& all, const int alpha, double& convergence)
{
    // until when to do bethe-takahashi steps only
    const double convergence_bt = 1e-3;
    // at what value of 'convergence' to stop and return an error
    const double divergence_break = 1e+4;
    // minimum number of steps to do without deviations
    const int min_iter_bt = 2;

    ++iterations_;

     if (!have_inner_sign_) {
        new_lambda_ = stepRapidity(all, alpha);

        if ( iterations_ > min_iter_bt && sq(lambda_-new_lambda_) < convergence_bt ) {
            have_inner_sign_ = true;
            setInitialDeviations(all, alpha);

            //set up for secant method
            last_lambda_ = lambda_;
            last_epsilon_ = epsilon_;
            last_delta_ = delta_;

            step_lambda_=new_lambda_;
            if (!stepDeviation(all, alpha, step_delta_, step_epsilon_))
                return IterResult::iter_err;

            // run a second step off step values
            lambda_ = step_lambda_;
            epsilon_ = step_epsilon_;
            delta_ = step_delta_;

            new_lambda_ = stepRapidity(all, alpha);
            if (!stepDeviation(all, alpha, new_delta_, new_epsilon_))
                return IterResult::iter_err;

            // set current values back
            lambda_ = last_lambda_;
            epsilon_ = last_epsilon_;
            delta_ = last_delta_;
        }

        convergence += sq(lambda_-new_lambda_);
        // don't signal convergence to calling function
        convergence+=1.;
        ++iterations_;
        return IterResult::iter_takahashi;
    }

    step_lambda_ = stepRapidity(all, alpha);

    // if not finite, nothing we can do.
    if (!finite(step_lambda_))
        return IterResult::iter_err_nonfinite;

    // simple step
    new_lambda_ = step_lambda_;
    convergence += sq(lambda_ - new_lambda_);

    if (convergence > divergence_break)
        return IterResult::iter_err_diverged;

    // 1-strings have no deviations
    if (string_length_==1)
        return IterResult::iter_ok;

    vector<double> new_step_delta (step_delta_.size());
    vector<double> new_step_epsilon (step_epsilon_.size());

    // find deviations & calculate their convergence. if not found, give up.
    if (!stepDeviation(all, alpha, new_step_delta, new_step_epsilon))
        return IterResult::iter_err;


    for (int a=0; a<string_length_;++a) {
        //new_epsilon_[a] = secantStep(last_epsilon_[a], epsilon_[a], step_epsilon_[a], new_step_epsilon[a]);
        //new_delta_[a] = secantStep(last_delta_[a], delta_[a], step_delta_[a], new_step_delta[a]);


        new_epsilon_[a] = dampedStep(epsilon_[a], new_step_epsilon[a], damping_dev_);
        new_delta_[a] = dampedStep(delta_[a], new_step_delta[a], damping_dev_);
        step_delta_ = new_step_delta;
        step_epsilon_ = new_step_epsilon;
    }

    // enforce holds
    for (int a=0; a<string_length_;++a) {
        if (hold_[a]) {
            new_delta_[a] = delta_[a];
            new_epsilon_[a] = epsilon_[a];
        }
    }

    // calculate convergence
    for (int a=0; a < string_length_; ++a)
        convergence += sq(new_delta_[a] - delta_[a]) + sq(new_epsilon_[a] - epsilon_[a]);

    if (convergence > divergence_break)
        return IterResult::iter_err_diverged;

    // report success
    return IterResult::iter_ok;
}



void NonsymString::refresh()
{
    last_lambda_ = lambda_;
    last_epsilon_ = epsilon_;
    last_delta_ = delta_;

    lambda_ = new_lambda_;
    epsilon_ = new_epsilon_;
    delta_ = new_delta_;
}


IterResult NonsymString::initiate(const vector<NonsymRoots*>& all, const int alpha)
{
    number_roots_ = 0;
	for (int beta=0; beta<all.size();++beta) {
        number_roots_ += all[beta]->size();
    }
    // if lambda not set, set non-interacting solution
    if (lambda_==0.) lambda_ = 0.5*string_length_*tan(0.5*ix2_*PI/big_n_);

    return IterResult::iter_ok;
}


void NonsymString::setInitialDeviations(const std::vector<NonsymRoots*>& all, const int alpha)
{

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
        inner_narrow_ = innerPairIsNarrow(all, alpha);
        if (inner_narrow_) {
            delta_[a_mid] = -initial_deviation_;
            delta_[a_mid+1] = initial_deviation_;
        }
        epsilon_[a_mid+1] = 0.;
    }
}
