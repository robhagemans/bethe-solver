#ifndef GAP_BASE_H
#define GAP_BASE_H

#include "base.h"


class Gap_Base: public Base {
public:
	Gap_Base (const Chain* const on_chain, const BaseData& base_data);
	Gap_Base (const Chain* const on_chain, const vector<int>& structure, const int number_holes=0);
	Gap_Base (const Chain* const on_chain, const int new_number_down, const vector<int>& structure, const int new_number_holes);
	// ground base
	Gap_Base (const Chain* const on_chain, const int the_number_down);
	// copy constructor
	Gap_Base (const Base& original): Base(original)
	{ 	if (original.p_chain->delta() <= 1.0) throw Exception("GapBase copy", exc_CopyMismatch); }
	// destructor
	virtual ~Gap_Base (void) {};
	// clone
	virtual Base* clone(void) const { return new Gap_Base(*this); };
	// ceil(2*I_max)   limit == maximum +1
	vector<int> limQuantumNumbers (void) const;
};


#endif
