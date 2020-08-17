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
	const REAL deviation_threshold = DEFAULT_DEVIATION_THRESHOLD,
	const long long int id_start=NO_ID, const long long int id_stop=NO_ID
);




/** 'scan' a single state **/
bool NEWscanState(
	AddFunc& addfunc,
	Quantity& quantity,
	State* p_state,
	const REAL deviation_threshold,
	const Policy policy,
	ScanResult& scan_result)
{

	State* p_right_state = 0;
	REAL deviation;

	try {
		// solve the Bethe Equations up to convergence
		if (p_state->solve(policy)) {
			// our state has converged
			deviation = p_state->stringDeviation();
			// our core business: calculate form factor
			if (deviation < deviation_threshold)
				quantity.setRightState(p_state);
			else {
 				p_right_state = new XXXDeviatedState ( *p_state );
 				if (p_right_state->solve(policy)) {
					quantity.setRightState(p_right_state);
				}
				else {
					delete p_right_state;
					return false;
				}
			}
		}
		else {
			addfunc.log() << name(p_state->p_base) <<"_state"<< p_state->id() <<": bethe equations not converged"<<endl;
			return false;
		}

		// calculate form factor and do whatever we want with it
		addfunc(quantity);

	}
	catch (Exception exc) {
		addfunc.log() << name(p_state->p_base) <<"_state"<< p_state->id() <<": "<< exc <<endl;
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
	State* p_state = newState (p_base, NO_ID);
	// loop over given interval of base
	for (long long int id = start_id; id < stop_id; ++id) {
		// set the id. may throw exceptions, which we don't catch as we can't resume anyway.
		int	state_change = p_state->setId(id);
		if (state_change>1 || !scan_result.number_total) p_state->setFreeRapidities();
		NEWscanState(addfunc, quantity, p_state, deviation_threshold, policy, scan_result);
	}

	// sum rule
	REAL max_sum = quantity.maxSum ();

	delete p_state;

	return scan_result.sum/max_sum;
}







class Counter : public OnlyLog {
public:
	long long int count;
	vector<long long int> deviate_count;
	REAL sum;
	inline Counter (void) : OnlyLog(cerr), count(0), sum(0.0), deviate_count(10) { } ;
	inline void operator() (Quantity& quantity)	{
		++count;
		REAL dev_mag = quantity.pRightState()->devianceMagnitude();
		int log_dev = 0;
		if ((dev_mag!=0)) {
			log_dev = 9 + (int) trunc(::log10(dev_mag));
			if (log_dev<0) log_dev=0;
			if (log_dev>9) log_dev=9;
		}
		sum += dev_mag;
		++deviate_count[log_dev];
	};
};





int run(void)
{
	int number_sites_max, number_sites_min;
	int number_down_left, number_down_diff;
	Policy policy = {2000, 1, 1.5* 1e-24, 1e-2, 1e-2, 10.0, 7, 20};

	cerr<<"N_min   N_max  N-M "<<endl;
	cin >> number_sites_min >>  number_sites_max >> number_down_diff;


	stringstream file_name;
	file_name << "N"<<number_sites_min<< "--"<<number_sites_max<<"_N-M"<<number_down_diff<<".count";
	const char* outfile_name = file_name.str().c_str();
	ofstream outfile (outfile_name);

	for (int number_sites=max(2*number_down_diff+2,number_sites_min); number_sites<number_sites_max; number_sites+=2)
	{
		number_down_left = number_sites/2 - number_down_diff;
		policy.precision = 1.5*1e-24*number_down_left;
		Chain* p_chain = newChain (1.0, number_sites, /* cutoff_types = */ 8);
		Base* p_ground_base = newGroundBase (p_chain, number_down_left);
		State* p_ground_state = newGroundState (p_ground_base);
		Quantity* p_quantity = newQuantity(0, p_ground_state);

		p_ground_state->solve();

		int number_down_right = p_quantity->rightNumberDown();
		vector< BaseData > all_bases = allNewBases(p_chain, number_down_right, /* max_string_length= */ number_down_right, /* max_number_particles= */ number_down_right, /*max_number_spinons= */ number_down_right, /*max_infinite= */ 0);//p_quantity->maxInfinite());
		Stopwatch calculation_time;

		Counter counter;

		REAL contrib=0.0;
		for (int i=0; i<all_bases.size(); ++i) {
			Base* p_base = newBase(p_chain, all_bases[i]);
			contrib += NEWscanBase (counter, *p_quantity, p_base, policy, /* deviation_threshold= */ 1e-20);
			delete p_base;
		}

		int expected_number = choose(number_sites, number_down_right) - choose(number_sites, number_down_right-1);

		outfile<<number_sites<<SEP<<number_down_left<<SEP;
		outfile<<expected_number<<SEP;
		outfile<<counter.count<<SEP;
		outfile<<-counter.count+ expected_number<<SEP;
		outfile<<counter.sum/counter.count<<SEP;
		outfile<<counter.deviate_count<<SEP;

		outfile<<endl;
	}
	outfile.close();
}
