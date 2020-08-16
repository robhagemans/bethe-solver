#include "chain.h"

const char* exc_OddChain = "number of sites must be even";
const char* exc_CopyMismatch = "argument mismatching in copy constructor";

/** copy constructor **/
Chain::Chain (const Chain& original) : anisotropy (original.anisotropy), number_sites (original.number_sites), number_types(original.number_types)
{
	string_length = new int[number_types];
	string_parity = new int[number_types];
	for (int i=0; i< number_types; ++i ) {
		string_length[i] = original.string_length[i];
		string_parity[i] = original.string_parity[i];
	}
};


/** destructor **/
Chain::~Chain(void)
{
	delete[] string_length;
	delete[] string_parity;
}

