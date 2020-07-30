#include "iso-base.h"

// set this flag if we want to use the formula from our paper instead of Takahashi's. It's including too many states, however.
// Takahashi was right after all, we just don't understand why.
//#define DO_NOT_USE_TAKAHASHIS_FORMULA

	/** construct from BaseData **/
	XXX_Base::XXX_Base (const Chain* const on_chain, const BaseData& base_data) 
		: Base(on_chain, base_data), number_infinite(0)
	{
		int number_roots = 0;
		for (int i=0; i< numberTypes(); ++i) number_roots += numberStringsOfType(i) * p_chain->stringLength(i);
		number_infinite = number_magnons - number_roots;

		calculateNumberSlots ();
		setNumberSpinons(number_spinons); 
		checkConsistency();
	}
	
	
	/** given-base constructor  **/
	XXX_Base::XXX_Base (const Chain* const on_chain, const vector<int>& structure, const int new_number_spinons, const int number_infinite_rapidities)
		: ::Base (on_chain, structure, new_number_spinons), number_infinite(number_infinite_rapidities)
	{
		calculateNumberSlots ();
		setNumberSpinons(new_number_spinons);
		checkConsistency();
		// the field number_roots does not contain the actual number of down spins! rather, it contains the number of finite roots
		//number_roots += number_infinite_rapidities; 
	}	
	
	/** constructor: ignores first element of structure, calculates from number_down and number_infinite  **/
	XXX_Base::XXX_Base (const Chain* const on_chain, const int new_number_down, const vector<int>& structure, const int new_number_spinons, const int number_infinite_rapidities)
		: ::Base (on_chain, (new_number_down-number_infinite_rapidities), structure, new_number_spinons), number_infinite(number_infinite_rapidities)
	{
		// NOTE: we inited by fooling parent into thinking the number of magnons was actually our number of roots. now we have to correct.
		number_magnons = new_number_down;
		
		calculateNumberSlots ();
		setNumberSpinons(new_number_spinons);
		checkConsistency();
	}
	
	/** ground-base constructor **/
	XXX_Base::XXX_Base (const Chain* const on_chain, const int the_number_roots)
		: ::Base(on_chain, the_number_roots), number_infinite(0)
	{
		calculateNumberSlots ();
		setNumberHoles(0);
		checkConsistency();
	}	
	

	/** calculate the boundary for allowed values of the half-integer quantum number I **/
	// for each string type j. lim_quantum == ceil(2*I_infty)
	// so that quantum_number < lim_quantum (strictly lower!)
	vector<int> XXX_Base::limQuantumNumbers (void) const
	{
		vector<int> lim_quantum(numberTypes());
		for (int j=0; j < numberTypes(); ++j) {
			int sum_theta = 0;
			//sum_theta += (numberStringsOfType(0)  - (j==0) ) * (1 + (j!=0) );
			for (int k=0; k < numberTypes(); ++k) 
#ifdef DO_NOT_USE_TAKAHASHIS_FORMULA
				// this sum_theta is as we thought it should be, based on 'lim x->infty (x-x) = 0'
				// equals (numberStringsOfType(k) - (j==k) ) * (2*min(chain.stringLength(j),chain.stringLength(k)) - (j==k) ), as in paper
				sum_theta += (numberStringsOfType(k) - (j==k) ) * (1 + 2*(min(p_chain->stringLength(j),p_chain->stringLength(k))-1) + (j!=k) );
#else		
				// this sum_theta is Takahashi's formula, which has the correct counting:
				sum_theta += numberStringsOfType(k) * (1 + 2*(min(p_chain->stringLength(j),p_chain->stringLength(k))-1) + (j!=k) );
#endif
			lim_quantum[j] = abs(p_chain->length() - sum_theta);
			// exclude the quantum numbers that yield 
			// infinite rapidities (if included, they would be double-counted anyway)
			// NOTE that quantum_number = 2*I and we have the limit *one* higher than the highest acceptable quantum_number (i.e. 1/2 higher in terms of I)
#ifdef DO_NOT_USE_TAKAHASHIS_FORMULA
			// NOTE also that if we use Takahashi's formula we already exclude them, so we should not use this.
 			--lim_quantum[j];
#endif
		}
		return lim_quantum;
	}

