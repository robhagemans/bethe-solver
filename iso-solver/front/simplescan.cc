
#include "simplescan.h"


using namespace std;



/** calculate the boundary for allowed values of the half-integer quantum number I **/
// for each string type j. lim_quantum == ceil(2*I_infty)
// so that quantum_number < lim_quantum (strictly lower!)
const int limQuantumNumber (const int big_n, const vector<int>& base, const int type)
{
	int sum_theta = 0;
	for (int k=0; k < base.size(); ++k)
		// this sum_theta is Takahashi's formula, which has the correct counting:
		sum_theta += base[k] * (1 + 2*(min(length(type),length(k))-1) + (type!=k) );

	return abs(big_n - sum_theta);
}



SimpleScanner::SimpleScanner(const int big_n, const int big_m)
    : big_n_(big_n), big_m_(big_m), base_(), ix2_()
{
}


const int SimpleScanner::nextNumber (const int number)
{
    if (number>0) return -number;
    else return -number+2;
}


bool SimpleScanner::firstInType (const int j, int max_alpha)
{
    if (max_alpha<0) max_alpha = base_[j];

    int next_qn = 1-base_[j]%2;
    for (int alpha=0; alpha < max_alpha; ++alpha) {
        // this happens only if it's a bad configuration - too many strings for a given type
        if (abs(next_qn) >= limQuantumNumber(big_n_, base_, j)) return false;

        ix2_[j][alpha] = next_qn;
        next_qn = nextNumber(next_qn);
    }
    return true;
}


bool SimpleScanner::firstState ()
{
 	ix2_.clear();
 	ix2_.resize(base_.size());

 	for (int j=0; j < base_.size(); ++j) {
        ix2_[j].resize(base_[j]);

        // return false if bad config (too many strings)
        if (!firstInType(j)) return false;
    }

    return true;
}


bool SimpleScanner::nextState ()
{
 	for (int j=0; j < base_.size(); ++j) {
 	    int lim_j = limQuantumNumber(big_n_, base_, j);

        // set ground for types less than j
        for (int k=0; k<j; ++k) {
            firstInType(k);
        }

        for (int alpha=0; alpha < base_[j]; ++alpha) {
            int next = nextNumber(ix2_[j][alpha]);
            if (abs(next) >= lim_j)
	            continue;
	        else if (alpha<base_[j]-1 && next == ix2_[j][alpha+1] )
	            continue;
	        else {
	            ix2_[j][alpha] = next;

	            // reset strings less than alpha
                firstInType(j, alpha);
                return true;
            }

        }
    }
    return false;
}




bool SimpleScanner::firstBase ()
{
    // all particles in first type
    base_ = { big_m_ };
    return true;
}


bool SimpleScanner::nextBase ()
{
    vector<int> old_base = base_;
    /*
    int number_down = 0;
    for (int j=0; j < base_.size(); ++j) {
        number_down += base_[j]*length(j);
    }
    */
    int number_down = big_m_;

    for (int j=1; j < number_down; ++j) {
        // try to increase the number of strings in type j, and set all other strings up to type j-1 to zero
        // use the real rapidities to fill up to the required number of particles

        if (j>=base_.size()) base_.resize(j+1);
        ++base_[j];
        base_[0] -= length(j); // string length

       // set all lower strings to zero
        for (int k=1; k < j; ++k) {
            base_[0] += length(k)*base_[k];
            base_[k] = 0;
        }

        // check if base is good & initialise the first state on the go
        if (base_[0]>=0 && firstState())  return true;

        // not good -> retain the old base and keep trying with next string type
        base_ = old_base;
    }
    return false;
}


