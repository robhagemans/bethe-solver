
#include "xxz-chain.h"

const char* exc_Rational = "rational zeta/pi";

//namespace Gapless {
	
	/** XXZ Chain: implementation **/
	
	
	XXZ_Chain::XXZ_Chain( const REAL delta_anisotropy, const int chain_length, const REAL precision, const int cutoff_fraction)
	{
		const char* here = "Chain::Chain";
		
		number_sites = chain_length;
		anisotropy = acos(delta_anisotropy);
		// my part_quot[i] == takahashi nu[i]
		// my denominator[i] == y[i] takahashi for i>=0
		// my sector == i takahashi
		// my part_sum[i] == m[i] takahashi
		
		if (number_sites%2) throw Exception (here, exc_OddChain); 
		
		REAL gammaoverpi = anisotropy/PI;
		
		vector<int> part_quot = findContinuedFraction (1.0/gammaoverpi, precision, cutoff_fraction); 
		// using 1/gammaoverpi == p0 inserts part_quot[0] == 0 (the integer part of the fraction is zero)
		// and shifts the rest of the fraction.
		// in this way my part_quot[i] == takahashi nu[i] even if I count from zero.

		// find denominators of continued fraction by Euler's recurrence relation
		// and find the partial sums of the sequence of partial quotients
		int number_quots = part_quot.size();

		// we cannot (yet) control rational zeta/pi. answers can and will be wrong. give up.
		if (number_quots < 3) throw Exception (here, exc_Rational);
		
		vector<int> denominator (number_quots); 
		vector<int> part_sum (number_quots); // sums of partial quotients
		denominator[0] = 1;
		part_sum[0] = 0;
		denominator[1] = part_sum[1] = part_quot[1];
		for (int i=2; i < number_quots; ++i) {
			denominator[i] = denominator[i-2] + part_quot[i] * denominator[i-1];
			part_sum[i] = part_sum[i-1] + part_quot[i];
		}
		// the number of allowed string types is the last of the partial sums,
		// aka the full sum, plus one. 
		number_types = part_sum.back()+1;
		// use the same cutoff as a max number of types.
		if (number_types > cutoff_fraction) number_types = cutoff_fraction;
		string_length = new int[number_types]; 
		string_parity = new int[number_types]; 
		
		for (int sector = 0; sector < number_quots-1; ++sector)
		{
			// my string_type == takahashi j, sort of, because my order differs. see notes.
			for (int string_type = part_sum[sector] + 1 ; string_type <= min (number_types, part_sum[sector+1]); ++string_type) {
				if (sector > 0) string_length[string_type] = denominator [sector-1]; // takahashi y[-1] == 0
				string_length[string_type] += (string_type - part_sum[sector])*denominator[sector];
				// parity defined by parity of floor((n_j-1)/p0)
				string_parity[string_type] = (int(floor((string_length[string_type]-1)*gammaoverpi)) % 2) ? -1 : 1;
			}
		}
		
			
		// after the loop, string_type 0 is undefined and 1 points to the positive-parity one string.
		// we'll redefine to add the negative-parity one-string as string_type 1,
		// which is always present and has to be added by hand anyway.
		string_length[0] = string_length[1] =1;
		string_parity[0] = 1;
		string_parity[1] = -1;
	}
	
	//XXZ_Chain::~XXZ_Chain(void)	{	}
	
	
		
//}
