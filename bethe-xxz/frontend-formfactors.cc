/* Rob Hagemans, instituut voor theoretische fysica, universiteit van Amsterdam */

#include <vector>
#include <iostream>
using namespace std;

#include "recursive.h"

// nullstream, thanks to Dietmar Kuehl, see http://www.velocityreviews.com/forums/t279756-null-output-stream.html
struct nullstream: std::ostream { nullstream(): std::ios(0), std::ostream(0) {} };
nullstream cnull;

#define CHOPHOLD 1e-10
inline complex<REAL> chop (complex<REAL> choppand)
{ 
	if (abs(real(choppand))< CHOPHOLD) real(choppand) = 0.0;
	if (abs(imag(choppand))< CHOPHOLD) imag(choppand) = 0.0;
}

// dump roots to the screen
class DumpRoots : public AddFunc {
public:
	inline DumpRoots(void) : AddFunc(cerr) { } ;
	inline void operator() (Quantity& quantity)	{ 
		State& state = quantity.rightState();
		cout << name(state.base) <<"_id"<< state.id() <<SEP<<quantity.formFactor()<<SEP<<quantity.energy()<<SEP<<state.norm()<<SEP<<
		state.quantum_number<<SEP<< state.roots()<<endl; 

		// ouch.
		XXX_State& iso_state = *( (XXX_State*) &state);
		vector<complex<REAL> > bethe_j = iso_state.calculateBetheI();
		REAL error = 0.0;
		for (int i=0; i<bethe_j.size();++i) {
			REAL two_j = 0.5*round(2.0*real(bethe_j[i]));
			error += norm(bethe_j[i]-two_j);
			cout<<two_j<<SEP;
		}
		cout<<endl;
		cout<<error<<endl;
		cout<<endl;
	};
};



int run(void)
{	
	double delta;
	int number_sites;
	int number_down_left, number_spinons;
	int number_infinite_rapidities;
	int function_index, element;
	int number_energy;
	REAL max_energy;
	
	//Policy policy = {1000, 1, 1e-28, 1e-2, 1e-2, 10.0, 1, 20}; // fewer newtons for large N; DEFAULT_POLICY is better for small systems
	Policy policy = {800, 100, 1e-30, 1e-2, 1e-2, 10.0, 7, 20};

	cerr<<"function  delta  N  M energy_bins max_energy"<<endl;
	cin >> function_index >> delta >> number_sites >> number_down_left >> number_energy >> max_energy;
	
	Chain* p_chain = newChain (delta, number_sites, /* cutoff_types = */ 8);	
	Base* p_ground_base = newGroundBase (*p_chain, number_down_left);
	State* p_ground_state = newGroundState (*p_ground_base);
	Quantity* p_quantity = newQuantity(function_index, *p_ground_state);
	p_ground_state->solve();

	vector< Base* > all_bases = allNewBases(p_chain, p_quantity->rightNumberDown(), /* max_string_length= */ 20, /* max_number_particles= */ 20, /*max_number_spinons= */ 20, p_quantity->maxInfinite());
	
	Stopwatch calculation_time;
// 	ScanIntervals interval;
	Matrix<REAL> basket (number_energy, number_sites);
	AddToMatrixFunc bin_this (basket, max_energy);
	AddToFileFunc dump_screen (cout, cerr); 
	DumpRoots dump_roots;
	cout.precision(30);
	
	REAL contrib=0.0;
	for (int i=0; i<all_bases.size(); ++i) {
 		contrib += scanBase (dump_roots, *p_quantity, all_bases[i], policy, acceptHalfBrillouin, /* deviation_threshold= */ 1e-20);
 		cerr<<"*** contribution "<<contrib<<endl;
	}	

}
