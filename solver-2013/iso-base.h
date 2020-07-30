#ifndef ISO_BASE_H
#define ISO_BASE_H

#include <vector>
using namespace std;

#include "base.h"

	class XXX_Base: public Base {
	public:
		int number_infinite;
		XXX_Base (const Chain* const on_chain, const BaseData& base_data);
		XXX_Base (const Chain* const on_chain, const vector<int>& structure, const int number_holes=0, const int number_infinite_rapidities=0);	
		XXX_Base (const Chain* const on_chain, const int new_number_down, const vector<int>& structure, const int new_number_holes, const int number_infinite_rapidities);
		// ground base
		XXX_Base (const Chain* const on_chain, const int the_number_roots);	
		// copy constructor
		XXX_Base (const Base& original): Base(original), number_infinite(original.numberInfiniteRapidities()) 
		{ 	if (original.p_chain->delta() != 1.0) throw Exception("XXXBase copy", exc_CopyMismatch); }
		// destructor
		virtual ~XXX_Base (void) {};
		// clone 
		virtual Base* clone(void) const { return new XXX_Base(*this); };	
		// ceil(2*I_max)   limit == maximum +1	
		virtual vector<int> limQuantumNumbers (void) const; 

		virtual inline int numberRoots (void) const { return number_magnons - number_infinite; };	
	};


#endif
