#include "xxz-base.h"

		
	XXZ_Base::XXZ_Base (const Chain* const on_chain, const BaseData& base_data) 
		: Base(on_chain, base_data)
	{
		const char* here = "XXZ_Base::XXZ_Base";
		
		int number_roots = 0;
		for (int i=0; i< numberTypes(); ++i) number_roots += numberStringsOfType(i) * p_chain->stringLength(i);
		if (number_magnons != number_roots) throw Exception(here, exc_Infinite);
	
		calculateNumberSlots ();
		setNumberSpinons(number_spinons); 
		checkConsistency();
	}

	/** given-base constructor  **/
	XXZ_Base::XXZ_Base (const Chain* const on_chain, const vector<int>& structure, const int new_number_spinons)
		: ::Base (on_chain, structure, new_number_spinons)
	{
		calculateNumberSlots ();
		setNumberSpinons(new_number_spinons);
		checkConsistency();
	}	
	
	XXZ_Base::XXZ_Base (const Chain* const on_chain, const int new_number_down, const vector<int>& structure, const int new_number_spinons)
		: ::Base (on_chain, new_number_down, structure, new_number_spinons)
	{
		calculateNumberSlots ();
		setNumberSpinons(new_number_spinons);
		checkConsistency();
	}
	
	
	/** ground-base constructor **/
	XXZ_Base::XXZ_Base (const Chain* const on_chain, const int the_number_roots)
		: ::Base(on_chain, the_number_roots)
	{
		calculateNumberSlots ();
		setNumberHoles(0);
		checkConsistency();
	}	
	
	//XXZ_Base::~XXZ_Base(void) {;}
	
	
	/** calculate the boundary for allowed values of the half-integer quantum number I **/
	// for each string type j. lim_quantum == ceil(2*I_infty)
	// so that quantum_number < lim_quantum (strictly lower!)
	vector<int> XXZ_Base::limQuantumNumbers (void) const
	{
		REAL zeta = p_chain->zeta();
		vector<int> lim_quantum(numberTypes());
		for (int j=0; j < numberTypes(); ++j) {
			REAL sum_theta = 0.0;
			int length_j = p_chain->stringLength(j);
			for (int k=0; k < numberTypes(); ++k) {
				REAL big_theta = 0.0;
				int length_k = p_chain->stringLength(k);
				if (p_chain->positiveParity(j)==p_chain->positiveParity(k)) {
					big_theta = 2.0 * atan( 1.0 / tan(  0.5*zeta*(length_j+length_k)  ) );
					if (length_j != length_k) 
						big_theta += 2.0 * atan( 1.0 / tan(  0.5*zeta*abs( length_j-length_k )) );
					for (int i=1; i< min(length_j, length_k); ++i) 
						big_theta += 2.0 * atan( 1.0 / tan(  0.5*zeta*(abs( length_j-length_k ) + 2*i)));
				}
				else {
					big_theta = -2.0 * atan( tan(  0.5*zeta*(length_j+length_k) ));
					if (length_j != length_k) 
						big_theta += -2.0 * atan( tan(  0.5*zeta*abs( length_j-length_k )  ));
					for (int i=1; i< min(length_j, length_k); ++i) 
						big_theta += -2.0 * atan( tan(  0.5*zeta*(abs(length_j-length_k ) + 2*i ) ));
				}
				sum_theta += numberStringsOfType(k) * big_theta; 
				if (j==k) sum_theta -= big_theta;
			}
			lim_quantum[j] = 
				int(ceil( (1.0/PI) 
					* fabs(
						p_chain->length()* 2.0 * 
						atan((	p_chain->positiveParity(j) 
								? (1.0/tan(  0.5*zeta*length_j)) 
								: (-tan(  0.5*zeta*length_j))
							))
						- sum_theta 
					) 
				));
		}
		return lim_quantum;
	}

