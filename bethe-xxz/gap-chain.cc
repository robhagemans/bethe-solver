
#include "gap-chain.h"


	Gap_Chain::Gap_Chain( const REAL delta_anisotropy, const int chain_length, const int cutoff_types)
	{
		const char* here = "Chain::Chain";

		number_sites = chain_length;
		anisotropy = acosh(delta_anisotropy);

		if (number_sites%2) throw Exception (here, exc_OddChain);

		number_types = cutoff_types;
		string_length = new int[number_types];
		string_parity = new int[number_types];
		for (int string_type= 0; string_type < number_types; ++string_type)	{
			string_length[string_type] = string_type+1;
			string_parity[string_type] = 1;
		}
	}

	Gap_Chain::~Gap_Chain(void)
	{
		delete[] string_length;
		delete[] string_parity;
	}

	inline REAL Gap_Chain::delta(void) const
	{ 	return cosh(anisotropy);	}
