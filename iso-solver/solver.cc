
#include "solver.h"

using std::vector;
using std::complex;
using std::min;
using std::cout;
using std::cerr;
using std::endl;
using std::abs;


bool isSymmetric(const vector< vector<int> >& ix2)
{
    for(int j=0; j<ix2.size(); ++j)
    for(int alpha=0; alpha<ix2[j].size(); ++alpha) {

        if (ix2[j][alpha]==0) continue;

        // find the opposite in the vector
        // we don't assume the vector is sorted so we have to check everything...
        bool found = false;
        for(int beta=0; beta<ix2[j].size(); ++beta) {
            if (alpha==beta) continue;
            if (ix2[j][beta] == -ix2[j][alpha]) {
                found = true;
                break;
            }
        }

        if (!found) return false;
    }
    return true;
}


// faster version where
// vector is assumed to be sorted 0 +1 -1 +2 -2 ...
bool isSymmetricSorted(const vector< vector<int> >& ix2)
{
    for(int j=0; j<ix2.size(); ++j)
    for(int alpha=0; alpha<ix2[j].size()-1; ++alpha) {
        // should happen only at alpha==0:
        if (ix2[j][alpha]==0) continue;
        // if there's an opposite, it must be the next element.
        if (ix2[j][alpha+1] != -ix2[j][alpha]) return false;
        // skip next element, it's our opposite
        ++alpha;
    }
    return true;
}


// guess sum of Bethe J's from Takahashi I configuration
// WARNING - this guess is often wrong!
// as it is, only used for symmetric configurations.
// to get Js, first solve the B-T equations to find the proper ordering - see class NonsymString
const int guessSum2xJ(const int big_n, const vector< vector<int> >& ix2, const int j, const int alpha)
{
	int l_j = length(j);  // we assume isotropic chain here
	int sum_2xj = ix2[j][alpha] + big_n*(l_j-1);
	for (int k=0; k < ix2.size(); ++k) {
		int l_k = length(k);
		for (int beta=0; beta < ix2[k].size(); ++beta) {
			if (j==k && beta==alpha) continue;

			// NOTE: here we assume the order of quantum numbers equals the ordering of rapidities
			// this is not necessarily true for Bethe Js - e.g. Kite state with J=-0.5 but positive lambda
			// is it true for Takahashi Is?
			// not sure, hence 'guess'
			// we use noninteracting rapidities: 0.5*length(j)*tan( 0.5*ix2[j][alpha]*PI/double(big_n)
			int sign = sgn( lambda0(big_n, ix2[j][alpha], j) - lambda0(big_n, ix2[k][beta], k) );

			// this works for central strings with 0 qns,
			// eseentially looking at the one with positive lambda
			if (!sign) sign = 1;

            sum_2xj -=  sign* ( l_j*l_k - 2*min(l_j, l_k) + (j==k) );

			//sum_2xj -= isgn(rapidity(j,alpha) - rapidity(k, beta))*(length_j*length_k - 2*min(length_j, length_k) + (j==k));
		}
	}
	return sum_2xj;
}


// string configuration constructor
// ix2 must be in the order 0, 1, -1, 3, -3, ...
StateSolver* newSolver(const int big_n, const std::vector<int>& base, const std::vector< std::vector<int> >& ix2, ReportFunc report)
{

    if (!isSymmetric(ix2)) {
        std::vector<NonsymRoots*> nonsymcomplexes;
        for(int j=0; j<ix2.size(); ++j)
        for(int alpha=0; alpha<ix2[j].size(); ++alpha)
            nonsymcomplexes.push_back(  new NonsymString( big_n, length(j), ix2[j][alpha] )  );
        return new Bunch<NonsymRoots>(nonsymcomplexes, report);
    }

    std::vector<SymRoots*> symcomplexes;
    std::vector<int> pos_real;
    std::vector<int> odd_zeros;
    std::vector<int> even_zeros;

    for(int alpha=0; alpha<ix2.at(0).size(); ++alpha) {
        if (ix2[0][alpha] == 0)
            odd_zeros.push_back(1);

        if (ix2[0][alpha] > 0)
            pos_real.push_back(guessSum2xJ(big_n, ix2, 0, alpha));
    }

    symcomplexes.push_back( new RealPairs(big_n,  pos_real ) );

    for(int j=1; j<ix2.size(); ++j)
    for(int alpha=0; alpha<ix2[j].size(); ++alpha) {

        if (ix2[j][alpha]==0) {
            if (length(j)%2) odd_zeros.push_back(length(j));
            else even_zeros.push_back(length(j));
        }
        else if (ix2[j][alpha]>0) {
            int jx2 = guessSum2xJ(big_n, ix2, j, alpha);
            symcomplexes.push_back( new SymString(big_n, length(j), jx2) );
        }
    }

    if (odd_zeros.size() == 0) ;
    else if (odd_zeros.size() == 1) {
        // only once central string or origin root
        symcomplexes.push_back( new CentralString(big_n, odd_zeros[0]) );
    }
    else if (odd_zeros.size()==2 && odd_zeros[0]==1 && odd_zeros[1]==3) {
        // 1&3 special kite
        symcomplexes.push_back( new Kite(big_n, 2) );
    }
    else if (odd_zeros.size() ==2) {

        // top and bottom pure imaginary roots
        std::vector<double> levels ((odd_zeros[1] - odd_zeros[0])/2);
        for (int i=0; i< levels.size(); ++i) {
            levels[i] = 1. + i + odd_zeros[0]/2;
        }
        symcomplexes.push_back( new ImagPairs(big_n, levels ));

        // get 2xJ for the shorter of the two central strings
        int jx2 = guessSum2xJ(big_n, ix2, odd_zeros[0]-1 , 0); /// WARNING - assumes length = type_nr+1

        symcomplexes.push_back( new SymString(big_n, odd_zeros[0], jx2) );

    }
    else if (odd_zeros.size() ==3) {
        symcomplexes.push_back( new CentralString(big_n, odd_zeros[0]) );

        // top and bottom pure imaginary roots
        std::vector<double> levels ((odd_zeros[2] - odd_zeros[1])/2);
        for (int i=0; i< levels.size(); ++i) {
            levels[i] = 1. + i + odd_zeros[1]/2;
        }
        symcomplexes.push_back( new ImagPairs(big_n, levels ));

        // get 2xJ for the shorter of the two central strings
        int jx2 = guessSum2xJ(big_n, ix2, odd_zeros[1]-1 , 0); /// WARNING - assumes length = type_nr+1

        symcomplexes.push_back( new SymString(big_n, odd_zeros[1], jx2) );

    }
    else{
        // construct generic kite state
        // not implemented
        return 0;
    }

    if (even_zeros.size() == 0) ;
    else if (even_zeros.size() == 1) {
        // only once central string or origin root
        symcomplexes.push_back( new CentralString(big_n, even_zeros[0]) );
    }
    else if (even_zeros.size() == 2) {

        // top and bottom pure imaginary roots
        std::vector<double> levels ((even_zeros[1] - even_zeros[0])/2);
        for (int i=0; i< levels.size(); ++i) {
            levels[i] = 0.5 + i + even_zeros[0]/2;
        }
        symcomplexes.push_back( new ImagPairs(big_n, levels ));
        // get 2xJ for the shorter of the two central strings
        int jx2 = guessSum2xJ(big_n, ix2, even_zeros[0]-1 , 0); // WARNING - assumes length = type_nr+1

        symcomplexes.push_back( new SymString(big_n, even_zeros[0], jx2) );

    }
    else if (even_zeros.size() == 3) {
        symcomplexes.push_back( new CentralString(big_n, even_zeros[0]) );

        // top and bottom pure imaginary roots
        std::vector<double> levels ((even_zeros[2] - even_zeros[1])/2);
        for (int i=0; i< levels.size(); ++i) {
            levels[i] = 0.5 + i + even_zeros[1]/2;
        }
        symcomplexes.push_back( new ImagPairs(big_n, levels ));
        // get 2xJ for the shorter of the two central strings
        int jx2 = guessSum2xJ(big_n, ix2, even_zeros[1]-1 , 0);  /// WARNING - assumes length = type_nr+1

        symcomplexes.push_back( new SymString(big_n, even_zeros[1], jx2) );

    }
    else {
        // construct generic kite state
        // not implemented
        return 0;
    }

    return new Bunch<SymRoots>(symcomplexes, report);

}





