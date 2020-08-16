#include "symmetric-string.h"

using namespace std;


/* CentralString */

CentralString::CentralString(const int big_n, const int string_length)
      : big_n_(big_n), string_length_(string_length), sum_jx2_(0),
        number_roots_odd_(false),
        delta_(string_length), new_delta_(string_length),
        hold_(string_length, false), knockout_(0), kite_mate_(0), string_mate_(0), mate_beta_(-1),
        iterations_(0)
{
    // use default policies, see header file
    setPolicy();
}

void CentralString::setPolicy(const double initial_deviation, const int max_simple_step, const double damping_delta)
{
    initial_deviation_ = initial_deviation;
    damping_delta_ = damping_delta;
    max_simple_step_ = max_simple_step;
}

bool CentralString::coupleKite(const Kite& mate, const int mate_beta)
{
    knockout_ = 3;//mate.string_length_;
    // set mate reference
    kite_mate_ = &mate;
    mate_beta_ = mate_beta;

    // mate must be shorter
    return (knockout_ <= string_length_);
}

bool CentralString::coupleString(const SymString& mate, const int mate_beta)
{
    knockout_ = mate.string_length_;
    // set mate reference
    string_mate_ = &mate;
    mate_beta_ = mate_beta;

    // mate must be shorter
    return (knockout_ <= string_length_);
}


vector<IVal> CentralString::getImagPairs() const
{
    // run through half string, excluding odd strings' centres
    const int n_pairs = (string_length_-knockout_)/2;
    vector<IVal> result (n_pairs);
    for (int a=0; a<n_pairs; ++a) {
        result[a] = { (string_length_-1) - 2*a, delta_[a] };
    }
    return result;
}

bool CentralString::hasOrigin() const
{
    return (!knockout_ && string_length_%2);
}



vector<CVal> CentralString::getRoots() const
{
    const int n_roots = string_length_ - knockout_;
    vector<CVal> result (n_roots);
    for (int a=0; a< (string_length_-knockout_)/2;  ++a) {
        result[a] =
            { 0.,0., (string_length_-1) - 2*a, delta_[a] };
        result[n_roots-1-a] =
            { 0.,0., -(string_length_-1) + 2*a, -delta_[a] };
    }
    if (!knockout_ && string_length_%2)
        result[string_length_/2] =  { 0.,0.,0,0. };
    return result;
}


int CentralString::size() const
{
    return string_length_ - knockout_;
}


IterResult CentralString::iterate(const vector<SymRoots*>& all, const int alpha, double& convergence)
{
	vector<double> step_delta;
	// secant method step for non-knocked out roots
    const int a_last = knockout_? (string_length_-knockout_)/2 : (string_length_-1)/2;  //exclude central pair even if no knockout
    const int a_mid = (string_length_-1)/2;

    // find deviations & calculate their convergence. if not found, give up.
    if (!stepCentralDeviation(all, alpha, step_delta))
        return IterResult::iter_err_neg_norm;

    for (int a=0; a<a_last; ++a)  {

        if (iterations_ < max_simple_step_)
            new_delta_[a] = dampedStep(delta_[a], step_delta[a], damping_delta_);
        else
            new_delta_[a] = secantStep(delta_last_[a], delta_[a], step_last_[a], step_delta[a]);

        new_delta_[string_length_-1-a] = - new_delta_[a];
    }

    // knocked out roots do simple step
    for (int a=a_last; a<=a_mid; ++a)  {
        new_delta_[a] = step_delta[a];
        new_delta_[string_length_-1-a] = - new_delta_[a];
    }

    // enforce holds
    for (int a=0; a<string_length_; ++a)
	    if (hold_[a]) new_delta_[a] = delta_[a];

    // calculate convergence
	for (int a=0; a < string_length_; ++a)
		convergence += sq(new_delta_[a] - delta_[a]);

    //delta_last_ = delta_;
    step_last_ = step_delta;
	//delta_ = new_delta;

    ++iterations_;
	return IterResult::iter_ok;
}

void CentralString::refresh()
{
    delta_last_ = delta_;
    delta_ = new_delta_;
    return;
}


void CentralString::setInitialValues(const double initial_deviation)
{
    for(int a=0; a<string_length_/2; ++a) {
        delta_[a] = initial_deviation * double(string_length_/2-a);
        delta_[string_length_-1-a] = -delta_[a];
    }
    // sets central root/inner pair to zero
    delta_[(string_length_-1)/2] = 0.;
    delta_[string_length_/2] = 0.;
}


IterResult CentralString::initiate(const vector<SymRoots*>& all, const int alpha)
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
    if (!stepCentralDeviation(all, alpha, step_last_))
        return IterResult::iter_err_neg_norm;
    delta_last_ = delta_;
    setInitialValues(2.*initial_deviation_);

    return IterResult::iter_ok;
}




/** calculate a new value for the deviance (delta) and aberration (epsilon) **/
bool CentralString::stepCentralDeviation (
                const vector<SymRoots*>& all, const int alpha,
                vector<double>& new_delta)
{
    // no deviations for origin root or pair.
	if (string_length_ < 3) 	{
	    new_delta.assign(string_length_, 0.);
		return true;
	}

	// l=5: 0 1 2 3 4, a_mid=2
    // l=4: 0 1 2 3, a_mid=1
    const int a_mid = (string_length_-1)/2;
    const int a_last = (string_length_-knockout_-1)/2;

    HPIVal special_theta = { 0, 0. };
	double special_rsq = 0.;

	double sum_log_rsq = 0.L;
	int sum_theta = 0; // number of PI terms
	vector<double> del_term (a_mid+1, 0.);
	double log_rsq_other = 0.L;

	for (int a=0; a <= a_last; ++a) {

        const IVal root = { (string_length_-1) - 2*a, delta_[a] };

        // scattering with others
        for (int beta=0; beta<all.size(); ++beta) {
            if (beta==alpha) continue;
	        if ((kite_mate_ || string_mate_) && beta==mate_beta_) continue;

            const vector<double> other_x = all[beta]->getRealPairs();
            const vector<IVal> other_y = all[beta]->getImagPairs();
            const vector<CVal> other_z = all[beta]->getComplexQuartets();

            int temp_theta =0;
            sum_log_rsq += 4.*sum_im_atan(root, all[beta]->hasOrigin(), other_x, other_y, other_z, temp_theta);
            sum_theta -= temp_theta/2;
        }

        // kinetic term
		//sum_theta += (abs(2.*y)>1.) ? big_n_*sgn(2.*y) : 0;
		//sum_theta += big_n_* re_atan_x4_div_pi( mul(root, 2) ) /2;
		sum_theta += big_n_* re_atan_x2_div_pi( mul(root, 2) );

		// 4.L*atih(2.L*y) == log( sq((y+0.5L) / (y-0.5L)) )
		sum_log_rsq -= big_n_* im_atan_x4( mul(root,2) );

        // contribution from 2*quantum number
		sum_theta -= 1-number_roots_odd_;

		// string self scattering minus dangerous terms
		for (int b=0; b < string_length_; ++b) {
		    // skip knocked out roots
		    if ( b>a_last && b<string_length_-1-a_last ) continue;
            // skip self
            if (b==a) continue;

            const IVal y_b = { (string_length_-1) - 2*b, delta_[b] };
            //? FIXME something is strange here. aren't we missing out half-terms?
	        //if (a-b!=+1 && a-b!=-1)  sum_theta -= re_atan_x4_div_pi( dif(root, y_b) )/2;
	        if (a-b!=+1 && a-b!=-1)  sum_theta -= re_atan_x2_div_pi( dif(root, y_b) );
	        //
	        if (a-b!=+1)  sum_log_rsq +=  zeta_plus_x2( dif(root, y_b) );
	        if (a-b!=-1)  sum_log_rsq -=  zeta_minus_x2( dif(root, y_b) );
		}


        if (kite_mate_) {
			// top kite root, excluding dangerous terms
	        const IVal y_top = { kite_mate_->pos_, kite_mate_->del_ }; //{ (string_length_-1) - 2*(a_last+1), kite_mate_->del_ };
	        //? FIXME: as above
	        //if (a!=a_last)  sum_theta -=  re_atan_x4_div_pi( dif(root, y_top) )/2;
		    if (a!=a_last)  sum_theta -=  re_atan_x2_div_pi( dif(root, y_top) );
		    //
		    sum_log_rsq +=  zeta_plus_x2( dif(root, y_top) );
		    if (a!=a_last)  sum_log_rsq -=  zeta_minus_x2( dif(root, y_top) );
	        // bottom kite root
	        //sum_theta -=  re_atan_x4_div_pi( sum(root, y_top) )/2;
		    sum_theta -=  re_atan_x2_div_pi( sum(root, y_top) );
		    sum_log_rsq += im_atan_x4( sum(root, y_top) );
		    // left & right kite roots
            sum_log_rsq +=  2.*im_atan_x4( sum(root, kite_mate_->x_) );
        }
        else if (string_mate_) {
            for (int k=0; k<string_mate_->string_length_; ++k) {

                const CVal z_k = { string_mate_->lambda_, string_mate_->epsilon_[k],
			                        (string_mate_->string_length_-1) - 2*k, string_mate_->delta_[k] };

			    // 2.* because both of   x+iy  -x+iy
			    // we include the (y+y_k) terms by running through the full string
                sum_log_rsq += 2.*zeta_plus_x2( dif(root, z_k) );

                // b > a and thus y>y_k as a is at most a_last
                // no theta effect as they're not on the imaginary axis
                if (a==a_last && 0==k) {
                    //special_rsq = -log( sq(z_k.real()) + sq(delta_[a] - string_mate_->delta_[k]) );
                    //special_theta += xi(z_k.real(), -delta_[a] + string_mate_->delta_[k]);
                    special_rsq = -zeta_minus_x2( dif(root, z_k) );
                    special_theta.add( xi_minus( dif(root, z_k) ) );
                }
                else  {
                    sum_log_rsq -= 2.*zeta_minus_x2( dif(root, z_k) );
                }
            }
        }

        // note: del_term[a] == new_delta[a] - new_delta[a+1]
        const double sign = (sum_theta%2) ? 1: -1;
	    del_term[a] = sign*exp(0.5L*sum_log_rsq);
	}

    new_delta.assign(string_length_, 0.);
    // inner pair has no deviations for central string

    bool success = true;

    if (kite_mate_)
        new_delta[a_last] = kite_mate_->del_ + del_term[a_last];

    if (string_mate_) {
        // positive sign choice
        del_term[a_last] = exp(0.5*sum_log_rsq)-sq(string_mate_->epsilon_[0] + string_mate_->lambda_);
        if (del_term[a_last] >= 0.) del_term[a_last] = sqrt(del_term[a_last]);
        else {
            // negative argument to sqrt
            success=false;
            // keep going?
        }

        new_delta[a_last] = string_mate_->delta_[0] + del_term[a_last];
        //new_delta[a_last] =-cos(0.5*sum_theta*PI+special_theta)*exp(0.5*(sum_log_rsq+special_rsq));
cerr<<" 1 "<<    string_mate_->delta_[0] + del_term[a_last] <<endl;
cerr<<" 2 "<<-cos(0.5*sum_theta*PI+special_theta.val())*exp(0.5*(sum_log_rsq+special_rsq)) <<endl;
cerr<<" 3 "<<-cos(0.5*sum_theta*PI+special_theta.val())*exp(0.25*sum_log_rsq) <<endl;
        //new_delta[a_last] = -cos(0.5*sum_theta*PI+special_theta)*exp(0.25*sum_log_rsq);
    }

    // add differences
    for (int a=a_last-1; a>=0; --a)
	    new_delta[a] += new_delta[a+1] + del_term[a];

	// set the conjugates
	// note that inner pair/string centre are necessarily zero -
    // have been zeroed on initialisation
    // new_delta[a_mid] = 0.;
	for (int a=a_mid+1; a<string_length_;++a)
		new_delta[a] = - new_delta[string_length_-a-1];

	return success;
}



