
#include "scan.h"


/** exceptions **/

const char* exc_InvalidIndex = "invalid momentum index";



/** policies **/

const int AddToFile::output_precision = 15;



/** addfunc **/

void AddToMatrixFunc::operator() (Quantity& quantity)
{
	const string here = "AddToMatrixFunc::operator()";
	int energy_resolution = table.height();
	
	int index_momentum = quantity.mode();
	if ((index_momentum < 0 ) || (index_momentum > table.width() )) throw Exception(here, exc_InvalidIndex);
	
	int index_energy = int(floor( 1.0 * energy_resolution * (quantity.energy() - min_energy) / (max_energy - min_energy) ));
	if ((index_energy > energy_resolution -1) || (index_energy < 0)) return;
	
 	table[index_energy][index_momentum] += quantity.formFactor();
}
	



/** 'scan' a single state **/
bool scanState(
	AddFunc& addfunc, 
	Quantity& quantity, 
	State* const p_state_in, 
	const REAL deviation_threshold,
	ScanResult& scan_result)
{
	State* p_state = p_state_in;
	// weights for half Brillouin zone
	// add to statistics only if the contribution is in the correct half-Brillouin sector
	// NOTE: this should give correct results for half-Brillouin as well as full-Brillouin, 
	// but may botch up for single momentum
	int weight = 0;
	if (p_state->mode()==0 || p_state->mode()==p_state->p_chain->length()/2) weight = 1;
	else if (p_state->mode()> 0  && p_state->mode() < p_state->p_chain->length()/2) weight = 2;
	// NOTE: hard-coded half Brillouin criterion:
	else return false;
	scan_result.number_total += weight;
	
	/// TODO: admissibility check: should be moved into State.solve(), ideally
	if (!p_state->admissible()) {
		scan_result.number_forbidden += weight;
#ifndef ADMIT_INADMISSIBLE
		addfunc.log(p_state)<<": inadmissible, mode "<<p_state->mode()<<" q.n. "<<p_state->quantum_number<<endl;
		return false;
#endif
	}

	try {
		quantity.setRightState(p_state);
		p_state = solve(p_state, deviation_threshold, quantity.deviable()); 
		quantity.setRightState(p_state);
	}
	catch (Exception exc) {
		if (exc.error == exc_NotConverged) scan_result.number_not_converged += weight;
		if (exc.error == exc_Deviated) scan_result.number_deviated += weight;
		quantity.setRightState(p_state_in);
		if (p_state != p_state_in) delete p_state;
		addfunc.log(p_state, exc);
		return false;
	}
	try {	
		// calculate form factor and do whatever we want with it
		addfunc(quantity);
		scan_result.number_calculated += weight;
		scan_result.sum += quantity.formFactor() * REAL(weight);
	} 
	catch (Exception exc) {
		if (exc.error == exc_Equal) scan_result.number_equal += weight;
		else scan_result.number_exceptions += weight;
		addfunc.log(p_state, exc);
		quantity.setRightState(p_state_in);
		if (p_state != p_state_in) delete p_state;
		return false;
	};
	return true;
}


		

/** scan through a base **/
REAL scanBase (
	AddFunc& addfunc, 
	Quantity& quantity, 
	Base* p_base, 
	const REAL deviation_threshold,
	long long int start_id, 
	long long int stop_id)
{
	const char* here = "scanBase";
	long long int lim_id = p_base->limId().back();
	if (start_id==NO_ID) start_id = 0;
	if (stop_id==NO_ID) stop_id = lim_id;
	

	// declare our calculation
	addfunc.log() <<"form_factor " << quantity.name() <<endl;
	addfunc.log() <<"delta "<< p_base->p_chain->delta()  <<endl;
	addfunc.log() <<"number_sites "<< p_base->p_chain->length() <<endl;
	addfunc.log() <<"left_number_down "<< quantity.leftNumberDown() <<endl;
	addfunc.log() <<"right_number_down "<< p_base->numberDown() <<endl;
	addfunc.log() <<"right_number_holes "<< p_base->numberHoles() <<endl;
	addfunc.log() <<"right_number_spinons "<< p_base->numberSpinons() <<endl;
	addfunc.log() <<"right_number_freedoms "<< p_base->numberFreedoms() <<endl;
	addfunc.log() <<"right_base "<< name(p_base) <<endl;
	addfunc.log() <<"0<=id< "<< lim_id << endl;
	
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
		scanState(addfunc, quantity, p_state, deviation_threshold, scan_result);
	}

	// log the scan result
	addfunc.log() <<"total "<<scan_result.number_total <<endl;
	addfunc.log() <<"calculated "<< scan_result.number_calculated <<endl;
	addfunc.log() <<"not_converged "<< scan_result.number_not_converged <<endl;
	addfunc.log() <<"deviated "<< scan_result.number_deviated <<endl;
	addfunc.log() <<"forbidden "<< scan_result.number_forbidden <<endl;
	addfunc.log() <<"equal "<< scan_result.number_equal <<endl;
	addfunc.log() <<"exceptions "<< scan_result.number_exceptions <<endl;
	
	// sum rule	
	REAL max_sum = quantity.maxSum ();
	
	// log the final result and time
	addfunc.log() <<"cpu_time " << stopwatch.humanReadable() <<endl;
	addfunc.log() <<"sum "<< scan_result.sum <<endl;
	addfunc.log() <<"maximum_sum "<< max_sum <<endl;
	addfunc.log() <<"fraction_sum "<< scan_result.sum/max_sum <<endl;
	addfunc.log(CODA);	
	
	delete p_state;
	return scan_result.sum/max_sum;
}

