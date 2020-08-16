#ifndef XXZ_BASE_H
#define XXZ_BASE_H

#include "base.h"

//namespace Gapless {

	class XXZ_Base: public Base {
	public:
		XXZ_Base (const Chain* const on_chain, const BaseData& base_data);
		// these constructors create a base with a given structure and number of holes
		XXZ_Base (const Chain* const on_chain, const vector<int>& structure, const int number_holes=0);
		// do same; calculate the first element of the base from the given number down
		// NOTE: caller must supply a dummy first base element in structure!
		XXZ_Base (const Chain* const on_chain, const int new_number_down, const vector<int>& structure, const int new_number_holes);
		// create a state with only particles in the 1+interval and no holes (ground state base)
		XXZ_Base (const Chain* const on_chain, const int the_number_down);
		// copy constructor
		XXZ_Base (const Base& original): Base(original)
		{ 	if (original.p_chain->delta() >= 1.0 || original.p_chain->delta() <= 0.0 ) throw Exception("XXZBase copy", exc_CopyMismatch); };
		// destructor
		virtual ~XXZ_Base (void) {};
		// clone
		virtual Base* clone(void) const { return new XXZ_Base(*this); };
		// ceil(2*I_max)   limit == maximum +1
		virtual vector<int> limQuantumNumbers (void) const;
	};
//}

#endif
