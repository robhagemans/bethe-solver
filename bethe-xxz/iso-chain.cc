#include "iso-chain.h"


	/** XXX Chain: implementation **/


	XXX_Chain::XXX_Chain(const int chain_length, const int cutoff_types)
	{
		number_sites = chain_length;
		number_types = min(chain_length, cutoff_types);
		if (number_sites%2) throw Exception ("XXX_Chain::XXX_Chain", exc_OddChain);

		// initialise superfluous fields
		anisotropy = 1.0;
		string_length = new int[number_types];
		string_parity = new int[number_types];
		for (int i=0; i < number_types; ++i) {
			string_length[i] = (i+1);
			string_parity[i] = 1;
		}
	}
