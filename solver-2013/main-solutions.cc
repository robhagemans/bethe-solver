/* Rob Hagemans, instituut voor theoretische fysica, universiteit van Amsterdam */

#include <vector>
#include <iostream>
using namespace std;

#include "generic.h"
#include "quantity.h"
#include "scan.h"

// dump roots to cout
class DumpRoots : public OnlyLog {
public:
	int count;
	REAL sum;
	inline DumpRoots() : OnlyLog(cerr), count(0), sum(0.0) { } ;
	virtual inline void operator() (Quantity& quantity)	{ 

		State* p_state = quantity.pRightState();

		cout<<p_state->quantum_number<<SEP;
		
		for (int j=0; j<p_state->p_base->numberTypes();++j)	
		for (int alpha=0; alpha< p_state->p_base->numberStringsOfType(j);++alpha) 
		for (int a=0; a<p_state->p_chain->stringLength(j); ++a)	{
			complex<REAL> ze_root = p_state->root(j, alpha, a);
			cout<<real(ze_root)<<SEP<<imag(ze_root)<<SEP;
		}

		cout<<"conv: "<<p_state->convergence<<SEP;
		
		XXX_State& iso_state = *( (XXX_State*) p_state);
 		
		cout<<"2J: "<<iso_state.calculateBetheI()<<SEP;
 		cout<<"error: "<<iso_state.betheError()<<SEP;
 		
		cout<<endl;
		++count;
	};
};

class DumpRootsForString: public DumpRoots {
public:
	int string_type_;
	inline DumpRootsForString(int string_type) : DumpRoots(), string_type_(string_type) { } ;
	virtual inline void operator() (Quantity& quantity)	{ 
	 
		State* p_state = quantity.pRightState();

		int j = string_type_; 
		cout<<p_state->quantum_number(j,0)<<SEP;
		cout<<p_state->quantum_number(j,0)%(2*(j+1))<<SEP;
		cout<<p_state->quantum_number<<SEP;
		
		for (int alpha=0; alpha< p_state->p_base->numberStringsOfType(j);++alpha) 
		for (int a=0; a<p_state->p_chain->stringLength(j); ++a)	{
			complex<REAL> ze_root = p_state->root(j, alpha, a);
			cout<<real(ze_root)<<SEP<<imag(ze_root)<<SEP;
		}

		cout<<"conv: "<<p_state->convergence<<SEP;
		
		XXX_State& iso_state = *( (XXX_State*) p_state);
 		
		cout<<"2J: "<<iso_state.calculateBethe2I()<<SEP;
 		cout<<"error: "<<iso_state.betheError()<<SEP;
 		
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

	//cerr<<"N  M"<<endl;
	//cin >> number_sites >> number_down_left;
	//number_sites = 22;
	//number_down_left = 2;
	number_sites = 48;
	number_down_left = 4;  //4
		
	Chain* p_chain = newChain (delta, number_sites, 8);	// cutoff_types=8
	Base* p_ground_base = newGroundBase (p_chain, number_down_left);
	State* p_ground_state = newGroundState (p_ground_base);
	Quantity* p_quantity = newQuantity(666, p_ground_state);
	
	// override standard policies and solve
	State::precision = 1e-26*number_sites*number_sites;
	
	State::max_iterations = 20000;
	p_ground_state->solve();

	int number_down_right = p_quantity->rightNumberDown();
	
	
	cout.precision(15);
	
#define ONESTATE
#define STRING_TYPE 3

//#define ALLBASE

//#define STRING_TYPE 1

#ifdef ALLBASE
// all bases:	
	DumpRoots dump_roots;
	vector< BaseData > all_bases = allNewBases(p_chain, number_down_right, /* max_string_length= */ number_down_right, /* max_number_particles= */ number_down_right, /*max_number_spinons= */ number_down_right, /*max_infinite= */ 0);//p_quantity->maxInfinite());
	for (int b=0; b< all_bases.size();++b) {
		Base* p_test_base = newBase(p_chain, all_bases[b]);
		scanBase (dump_roots, *p_quantity, p_test_base,  1e-20); // deviation_threshold= 1e-20	
	}
	
#else

// one base:		
	DumpRootsForString dump_roots(STRING_TYPE);
	// one n-string, all others real
	vector<int> base_vec (STRING_TYPE+1, 0);
	base_vec[0] = number_down_right-STRING_TYPE-1;
	base_vec[STRING_TYPE] = 1;
	
	// no spinons, just the string. 
	Base* p_test_base = newBase(p_chain, base_vec, 0);

#ifndef ONESTATE

  	scanBase (dump_roots, *p_quantity, p_test_base,  1e-20); // deviation_threshold= 1e-20

#else
// one state:
  	State* p_test_state = newState(p_test_base, 37);
  	//State* p_test_state = newState(p_test_base, 39);
  	ScanResult scan_result; 
  	scanState(dump_roots, *p_quantity, p_test_state, 1e-20, scan_result);

#endif
#endif

  	return 0;
}
