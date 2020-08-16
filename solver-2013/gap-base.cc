#include "gap-base.h"

	Gap_Base::Gap_Base (const Chain* const on_chain, const BaseData& base_data)
		: Base(on_chain, base_data)
	{
		const char* here = "Gap_Base::Gap_Base";

		int number_roots = 0;
		for (int i=0; i< numberTypes(); ++i) number_roots += numberStringsOfType(i) * p_chain->stringLength(i);
		if (number_magnons != number_roots) throw Exception(here, exc_Infinite);

		calculateNumberSlots ();
		setNumberSpinons(number_spinons);
		checkConsistency();
	}

	/** given-base constructor  **/
	Gap_Base::Gap_Base (const Chain* const on_chain, const vector<int>& structure, const int new_number_spinons)
		: ::Base (on_chain, structure, new_number_spinons)
	{
		calculateNumberSlots ();
		setNumberSpinons(new_number_spinons);
		checkConsistency();
	}


	Gap_Base::Gap_Base (const Chain* const on_chain, const int new_number_roots, const vector<int>& structure, const int new_number_spinons)
		: ::Base (on_chain, new_number_roots, structure, new_number_spinons)
	{
		calculateNumberSlots ();
		setNumberSpinons(new_number_spinons);
		checkConsistency();
	}

	/** ground-base constructor **/
	Gap_Base::Gap_Base (const Chain* const on_chain, const int the_number_roots)
		: ::Base(on_chain, the_number_roots)
	{
		calculateNumberSlots ();
		setNumberHoles(0);
		checkConsistency();
	}


	/** calculate the boundary for allowed values of the half-integer quantum number I **/
	// for each string type j. lim_quantum == ceil(2*I_infty)
	// so that quantum_number < lim_quantum (strictly lower!)
	vector<int> Gap_Base::limQuantumNumbers (void) const
	{
		vector<int> lim_quantum(numberTypes());
	/*
		for (int j=0; j < numberTypes(); ++j) {
			int sum_theta = 0;
			for (int k=0; k < numberTypes(); ++k)
				sum_theta += (numberStringsOfType(k) - (j==k) ) * (1 + 2*(min(chain.stringLength(j), chain.stringLength(k))-1) + (j!=k) );
			lim_quantum[j] = abs(chain.length() - sum_theta);
			// include the quantum numbers on the boundary
			// NOTE: this induces double counting; remove one of the boundary quantum numbers later!
			++lim_quantum[j];
		}
	*/
		// the very maximum of possibilities
		// many states will have to be scrapped as doubly-counted.
		for (int j=0; j < numberTypes(); ++j) lim_quantum[j] = p_chain->length()+1;
		return lim_quantum;
	}

