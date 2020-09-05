#ifndef BTSTRING_H
#define BTSTRING_H

#include <vector>
#include <complex>
#include <cmath>

#include "basics.h"


// string length for isotropic chain
inline const int length(const int type) { return type+1; }


// noninteracting rapidities
inline const double lambda0(const int big_n, const int ix2, const int j)
{
    return 0.5*length(j)*tan( 0.5*ix2*PI/double(big_n) );
}


/*
// extract sum of Bethe J's from Takahashi I configuration
const int sum2xJ(const int big_n, const int ix2, const double lambda_j, const int l_j, const vector<double>& lambdas, const vector<double>& lengths)
{
	int sum_2xj = ix2 + big_n*(l_j-1);
	for (int k=0; k < lambdas.size(); ++k) {
		const int l_k = lengths[k];
        sum_2xj -= isgn(lambda_j - lambdas[k]) * ( l_j*l_k - 2*min(l_j, l_k) + (l_j==l_k) );
	}
	return sum_2xj;
}
*/

#endif

