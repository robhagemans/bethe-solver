/* Rob Hagemans, instituut voor theoretische fysica, universiteit van Amsterdam */

#include <vector>
#include <iostream>
using namespace std;

#include "generic.h"
#include "quantity.h"
#include "scan.h"



/** scan through a full base and calculate form factors **/
REAL NEWscanBase (
	AddFunc& addFunc,
	Quantity& quantity,
	Base* p_base, 
	Policy policy, 
	const AcceptFunc& accept = acceptHalfBrillouin,
	const REAL deviation_threshold = DEFAULT_DEVIATION_THRESHOLD,
	const long long int id_start=NO_ID, const long long int id_stop=NO_ID
);




/** 'scan' a single state **/
bool NEWscanState(
	AddFunc& addfunc, 
	Quantity& quantity, 
	State* p_state, 
	const AcceptFunc& accept, 
	const REAL deviation_threshold,
	const Policy policy,
	ScanResult& scan_result)
{
	// weights for half Brillouin zone
	// add to statistics only if the contribution is in the correct half-Brillouin sector
	// NOTE: this should give correct results for half-Brillouin as well as full-Brillouin, 
	// but may botch up for single momentum
	int weight = 0;
 	if (p_state->mode()==0 || p_state->mode()==p_state->chain.length()/2) weight = 1;
 	else if (p_state->mode()> 0  && p_state->mode() < p_state->chain.length()/2) weight = 2;
	
	scan_result.number_total += weight;
	if (!accept(*p_state)) {
		scan_result.number_not_accepted += weight;
		addfunc.logstream<<name(p_state->base)<<"_id"<< p_state->id()<<": not accepted"<<endl;
		return false;
	}

	if (!p_state->admissible()) {
 		scan_result.number_forbidden += weight;
//  		addfunc.logstream<<name(p_state->base)<<"_id"<< p_state->id()<<": inadmissible, mode "<<p_state->mode()<<" q.n. "<<p_state->quantum_number<<endl;
//  		return false;
	}

	State* p_right_state = 0;
	REAL deviation;
					
								
	try {
		// solve the Bethe Equations up to convergence
		if (p_state->solve(policy)) { 
			// our state has converged
			deviation = p_state->stringDeviation();
			// our core business: calculate form factor
			if (deviation < deviation_threshold) quantity.setRightState(*p_state); 
			else if (
				p_state->chain.delta()==1.0 
				&& (
					(quantity.name() == Longitudinal::ff_name && p_state->base.numberInfiniteRapidities() == 0) 
					|| (quantity.name() == SPlusMinus::ff_name && p_state->base.numberInfiniteRapidities() < 2) 
				)
			) {
 				p_right_state = new XXXDeviatedState ( *p_state );
				if (!p_right_state->solve(policy)) {
					// we haven't converged! (counts as deviated)
					scan_result.number_deviated += weight;
					
					// log as exception (for debugging purposes)
					addfunc.logstream 
						<<name(p_right_state->base)<<"_id"<< p_right_state->id() 
						<<": deviance not converged after "<< p_right_state->iterations
						<<" iterations; q.n. "<<p_right_state->quantum_number<<endl;

					// clean up
					delete p_right_state;
					// I guess, God bless and goodbye
					return false;
				}
				else quantity.setRightState(*p_right_state);
			}
			
			else {
				scan_result.number_deviated += weight;
				return false;
			}
		}
		else {
			// we haven't converged!
			scan_result.number_not_converged += weight;
			return false;
		}
		// calculate form factor and do whatever we want with it
		addfunc(quantity);

		scan_result.number_calculated += weight;
// 		if (finite(quantity.formFactor()))	scan_result.sum += quantity.formFactor() * REAL(weight);
	} 
	catch (Exception exc) {
		if (exc.error == exc_Equal) scan_result.number_equal += weight;
		else scan_result.number_exceptions += weight;
		addfunc.logstream << name(p_state->base) <<"_state"<< p_state->id() <<": "<< exc <<endl;
	};
	
	// if we've created something new, delete it now.
	if (p_right_state && p_right_state != p_state) delete p_right_state;
	return true;
}


		






/** scan through a base **/
REAL NEWscanBase (
	AddFunc& addfunc, 
	Quantity& quantity, 
	Base* p_base, 
	Policy policy,  
	const AcceptFunc& accept, 
	const REAL deviation_threshold,
	long long int start_id, 
	long long int stop_id)
{
	const char* here = "scanBase";
	long long int lim_id = p_base->limId().back();
	if (start_id==NO_ID) start_id = 0;
	if (stop_id==NO_ID) stop_id = lim_id;
	
	// set timer
	Stopwatch stopwatch;
	// set sum to zero
	ScanResult scan_result; 
	// NO_ID: do not set id (and thus don't throw exceptions)
	State* p_state = newState (*p_base, NO_ID);	
	// loop over given interval of base
	for (long long int id = start_id; id < stop_id; ++id) {
		// set the id. may throw exceptions, which we don't catch as we can't resume anyway.
		int	state_change = p_state->setId(id); 
		if (state_change>1 || !scan_result.number_total) p_state->setFreeRapidities();	
		NEWscanState(addfunc, quantity, p_state, accept, deviation_threshold, policy, scan_result);
	}
	
	// sum rule	
	REAL max_sum = quantity.maxSum ();
	
	delete p_state;
// cerr<<"sums:     "<<scan_result.sum<<SEP<<max_sum<<SEP<<scan_result.sum/max_sum<<endl;	
	return scan_result.sum/max_sum;
}


























// nullstream, thanks to Dietmar Kuehl, see http://www.velocityreviews.com/forums/t279756-null-output-stream.html
struct nullstream: std::ostream { nullstream(): std::ios(0), std::ostream(0) {} };
nullstream cnull;

#define CHOPHOLD 1e-10
inline complex<REAL> chop (complex<REAL> choppand)
{ 
	if (abs(real(choppand))< CHOPHOLD) real(choppand) = 0.0;
	if (abs(imag(choppand))< CHOPHOLD) imag(choppand) = 0.0;
}

// dump roots to cout
class DumpRoots : public AddFunc {
public:
	int count;
	REAL sum;
	inline DumpRoots(void) : AddFunc(cerr), count(0), sum(0.0) { } ;
	inline void operator() (Quantity& quantity)	{ 
		State& state = quantity.rightState();
		cout<<count+1<<". &"<<SEP;
		
		// ouch.
		XXX_State& iso_state = *( (XXX_State*) &state);

		vector<complex<REAL> > bethe_j = iso_state.calculateBetheI();
 		REAL error = 0.0;
		for (int i=0; i<bethe_j.size();++i) {
			REAL two_j = 0.5*round(2.0*real(bethe_j[i]));
 			error += norm(bethe_j[i]-two_j);
  			cout<<"$"<<int(2.0*two_j)<<"$ &"<<SEP;
		}

 		
//   		cout<<sqrt(error)<<SEP;
// 		cout<<iso_state.betheError()<<SEP;
				
		for (int j=0; j< state.base.numberTypes();++j)
		for (int alpha=0; alpha< state.base.numberStringsOfType(j);++alpha) {
			int length_j = state.chain.stringLength(j);
			if (length_j>1){
				if ((j==state.base.numberTypes()-1) && alpha==state.base.numberStringsOfType(j)-1)
					cout<<"\\multicolumn{"<<length_j<<"}{c|}{$";
				else
					cout<<"\\multicolumn{"<<length_j<<"}{c}{$";
				cout<<state.quantum_number(j, alpha)<<"_"<<length_j<<"$} &"<<SEP;
			}
			else 
				cout<<"$"<<state.quantum_number(j, alpha)<<"_"<<length_j<<"$ &"<<SEP;
		}
		
		
		for (int j=0; j< state.base.numberTypes();++j)
		for (int alpha=0; alpha< state.base.numberStringsOfType(j);++alpha) 
		for (int a=0; a<(state.chain.stringLength(j)+1)/2; ++a)	{
			complex<REAL> ze_root = state.root(j, alpha, a);
			if (abs(imag(ze_root))<CHOPHOLD) {
				if (abs(real(ze_root))<CHOPHOLD) cout<<"$"<<0.0<<"$ &"<<SEP;
				else cout<<"$"<<real(ze_root)<<"$ &"<<SEP;
			}
			else {
				if ((j==state.base.numberTypes()-1) && alpha==state.base.numberStringsOfType(j)-1 && a==(state.chain.stringLength(j)+1)/2-1)
					cout<<"\\multicolumn{2}{c|}{$";
				else
					cout<<"\\multicolumn{2}{c}{$";
				if (abs(real(ze_root))<CHOPHOLD) cout<<"\\pm "<<abs(imag(ze_root))<<"i"<<"$} &"<<SEP;
				else cout<<real(ze_root)<<"\\pm "<<abs(imag(ze_root))<<"i"<<"$} &"<<SEP;
			}
			
		}
		
 		cout<<"$"<<state.energy()<<"$ \\\\"<<endl;
		
/*		
		XXXDeviatedState& dev_state = *( (XXXDeviatedState*) &state);
		try {
			if (dev_state.oddSymmetric()) {
				cerr<<name(state.base)<<"_id"<<state.id()<<SEP<<iso_state.betheError()<<SEP;
				cerr<<state.admissible()<<SEP<<dev_state.oddSymmetric()<<SEP<<state.symmetric()<<SEP;
				cerr<<quantity.formFactor()<<endl; 	
				cerr<<state.roots()<<endl<<state.quantum_number<<endl;;
				cerr<<dev_state.betheError()<<endl;
				
			}
 			if (finite(quantity.formFactor())) sum += quantity.formFactor();
		}
		catch (Exception exc) {	
 			cerr<<"(error)"<<endl;
		}
*/
		++count;
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
	Policy policy = {2000, 100, 1e-30, 1e-2, 1e-2, 10.0, 7, 20};

	cerr<<"delta  N  M"<<endl;
	cin >> delta >> number_sites >> number_down_left;
	
	Chain* p_chain = newChain (delta, number_sites, /* cutoff_types = */ 8);	
	Base* p_ground_base = newGroundBase (*p_chain, number_down_left);
	State* p_ground_state = newGroundState (*p_ground_base);
	Quantity* p_quantity = newQuantity(-1, *p_ground_state);
// 	Quantity* p_minus = newQuantity(1, *p_ground_state);
// 	Quantity* p_plus = newQuantity(-1, *p_ground_state);
	
	
	p_ground_state->solve();
/*
	//_R4_s0_base-1-0-1_id0:
	Base* p_base = newBase("_R4_s0_base-1-0-1", p_chain);
	State* p_state = newState(*p_base, 0);
	solve(p_state, policy,  1e-20);
	cout<<p_state->roots()<<endl;
	cout<<( (XXXDeviatedState*) p_state)->calculateBetheI()<<endl;
	return 0;
*/	
	int number_down_right = p_quantity->rightNumberDown();
	vector< Base* > all_bases = allNewBases(p_chain, number_down_right, /* max_string_length= */ number_down_right, /* max_number_particles= */ number_down_right, /*max_number_spinons= */ number_down_right, /*max_infinite= */ 0);//p_quantity->maxInfinite());
	
	Stopwatch calculation_time;
	
// 	Matrix<REAL> basket (number_energy, number_sites);
// 	AddToMatrixFunc bin_this (basket, max_energy);
// 	AddToFileFunc dump_screen (cout, cerr); 
	DumpRoots dump_roots;
	cout.precision(15);
	
	REAL contrib=0.0;
	for (int i=0; i<all_bases.size(); ++i) {
 		contrib += NEWscanBase (dump_roots, *p_quantity, all_bases[i], policy, acceptAlways, /* deviation_threshold= */ 1e-20);
	}	
	
	cerr<<endl;
	cerr<<contrib<<endl;
	cerr<<dump_roots.sum<<endl;
	int expected_number = choose(number_sites, number_down_right) - choose(number_sites, number_down_right-1);
	cout<<endl;
	cerr<<endl;
	cerr<<dump_roots.count<<SEP<<expected_number<<endl;
}
