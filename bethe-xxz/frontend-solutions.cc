/* Rob Hagemans, instituut voor theoretische fysica, universiteit van Amsterdam */

#include <vector>
#include <iostream>
using namespace std;

#include "generic.h"
#include "quantity.h"
#include "scan.h"

#define STRING_TYPE 3

// dump roots to cout
class DumpRoots : public OnlyLog {
public:
	int count;
	REAL sum;
	inline DumpRoots(void) : OnlyLog(cerr), count(0), sum(0.0) { } ;
	virtual inline void operator() (Quantity& quantity)	{ 
		State* p_state = quantity.pRightState();
		
		int j = STRING_TYPE; 
		
		cout<<p_state->quantum_number(j,0)<<SEP;
		cout<<p_state->quantum_number(j,0)%(2*(STRING_TYPE+1))<<SEP;
		
		for (int alpha=0; alpha< p_state->p_base->numberStringsOfType(j);++alpha) 
		for (int a=0; a<p_state->p_chain->stringLength(j); ++a)	{
			complex<REAL> ze_root = p_state->root(j, alpha, a);
			cout<<real(ze_root)<<SEP<<imag(ze_root)<<SEP;
		}
		cout<<p_state->convergence<<SEP;
		
		XXX_State& iso_state = *( (XXX_State*) p_state);
 		
		cout<<iso_state.calculateBethe2I()<<SEP;
 		cout<<iso_state.betheError()<<SEP;
 		
		cout<<endl;
		++count;
	};
};







int run(void)
{	
	double delta = 1.0;
	int number_sites;
	int number_down_left, number_spinons;
	int number_infinite_rapidities;
	int function_index, element;
	int number_energy;
	REAL max_energy;

/* calculate solutions for a single STRING_TYPE-string, filled up by real roots */

	cerr<<"N  M"<<endl;
	cin >> number_sites >> number_down_left;
	
	Chain* p_chain = newChain (delta, number_sites, 8);	// cutoff_types=8
	Base* p_ground_base = newGroundBase (p_chain, number_down_left);
	State* p_ground_state = newGroundState (p_ground_base);
	Quantity* p_quantity = newQuantity(666, p_ground_state);
	
	// override standard policies and solve
	State::precision = 1e-26*number_sites*number_sites;
	State::max_iterations = 2000;
	p_ground_state->solve();

	int number_down_right = p_quantity->rightNumberDown();
	
	// one n-string, all others real
	vector<int> base_vec (STRING_TYPE+1, 0);
	base_vec[0] = number_down_right-STRING_TYPE-1;
	base_vec[STRING_TYPE] = 1;
	
	// no spinons, just the string. 
	Base* p_test_base = newBase(p_chain, base_vec, 0);
	
	DumpRoots dump_roots;
	cout.precision(15);
	
  	scanBase (dump_roots, *p_quantity, p_test_base,  1e-20); // deviation_threshold= 1e-20
	
/* keep lambda constant and increase N */
/*
	int length_start, length_stop;

	cerr << "N_start N_stop"<<endl;
	cin >> length_start >> length_stop;
	
	REAL factor = 0.94;
	
	for (int number_sites = length_start; number_sites < length_stop; number_sites+=2) {
		
		int number_down_left = STRING_TYPE+1;
		
		Chain* p_chain = newChain (delta, number_sites, 8);	// cutoff_types=8
		Base* p_ground_base = newGroundBase (p_chain, number_down_left);
		State* p_ground_state = newGroundState (p_ground_base);
		Quantity* p_quantity = newQuantity(666, p_ground_state);
		
		// override standard policies and solve
		State::precision = 1e-26*number_sites*number_sites;
		State::max_iterations = 2000;
		p_ground_state->solve();
		
		
		int number_down_right = p_quantity->rightNumberDown();
		
		// one n-string, all others real
		vector<int> base_vec (STRING_TYPE+1, 0);
		base_vec[0] = number_down_right-STRING_TYPE-1;
		base_vec[STRING_TYPE] = 1;
		
		// no spinons, just the string. 
		Base* p_test_base = newBase(p_chain, base_vec, 0);
		State* p_test_state = newState(p_test_base, NO_ID);
		p_test_state->quantum_number(STRING_TYPE,0) = (int) round(factor*number_sites);
		
		
		try {
			p_quantity->setRightState(p_test_state);
			p_test_state = solve(p_test_state, 1e-20, true); 
			p_quantity->setRightState(p_test_state);
			
			DumpRoots dump_roots;
			cout.precision(15);
			dump_roots(*p_quantity);
		}
		catch (Exception exc) {
			cerr<<exc<<endl;
		}
		
	
	}
 /* */
}
