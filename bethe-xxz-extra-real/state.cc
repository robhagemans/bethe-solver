#include "bethe.h"
#include "state.h"
#include <stdio.h>

/** exceptions **/

const char* exc_NotGroundBase = "not a ground base ";
const char* exc_NegativeHoles = "negative number of holes";
//const char* exc_Forbidden = "forbidden state"; //State::SetID
const char* exc_Strings = "argument state must not have strings";
const char* exc_NumberDown = "wrong number of down spins in argument state";
const char* exc_Runaway = "runaway rapidity";
const char* exc_NonFinite = "non-finite result";
const char* exc_Uncontrol = "uncontrollable rapidity";
const char* exc_Equal = "equal rapidities on left and right of matrix element";	// matrixHm2P(), matrixHminus()
const char* exc_BadShift = "incompatible shifts";
const char* exc_NotImplemented = "not implemented";

/** sentinel values **/

const REAL NO_CONVERGENCE = sq(BOOMBANOOMBA)*4.0; // some big number.
const REAL NOT_CALCULATED = -NO_CONVERGENCE;


/** initialise statics of State **/

// default solution policy (works best for short chains)
// general
int State::max_iterations			= 1000;		// maximum iterations
REAL State::precision				= 1e-28;	// required precision for convergence
// newton's method
int State::max_newton 				= 100;		// maximum Newton's method steps
REAL State::newton_threshold		= 1e-2;		// precision to extrapolate to before newton's method is invoked
REAL State::newton_factor			= 1e-2;		// factor with which to increase newton threshold if newton doesn't work
REAL State::newton_bandwidth		= 10.0;		// factor of worsening we allow from newton's method before stopping it
// solve() strategy
int State::newton_consecutive		= 7;		// number of consecutive newton steps allowed
int State::extrapolate_consecutive  = 20;		// number of consecutive iterations/extrapolations allowed	
// policies for iterate() 
int State::run_max 					= 50;		// maximum number of steps when running
float State::run_sloth 				= 0.5; 		// between 0 and 1; higher number means slower run.


/** internal **/

// defines sorting order for quantum numbers (for use with std::sort() )
// true if a is less than b according to the ordering:
// 0, -2, 2, -4, 4, -6, 6, etc.
// -1, 1, -3, 3, -5, 5, etc.
inline bool quantumNumberLessThan (const int a, const int b) 
{ 	return (abs(a)<abs(b)) || ( (abs(a)==abs(b)) && (a<b) ); 	}

// fill quantum_number with the quantum numbers of the ground state
void setGroundQuantumNumbers(const Base* p_base, Strip<int>& quantum_number)
{
	// - boundary_0 is the largest integer lower than the first *quantum number* (2*I) in the ground sector.
	// (therefore, its parity is different from that of quantum numbers)
	// base should ensure it has the same parity as the number of slots it fits into
	int boundary_0 = p_base->number_slots_sector[0];
	// fill sector[0]. all holes on the boundary. 
	for (int instance = 1-p_base->number_slots_sector[0]%2; instance < p_base->number_particles_sector[0]; instance+=2) 
		quantum_number(0, instance) = instance;	
	for (int instance = p_base->number_slots_sector[0]%2; instance < p_base->number_particles_sector[0]; instance+=2) 
		quantum_number(0, instance) = - instance-1;
	// fill sector[1]. it consists of two parts.
	// it always has an even number of slots.
	for (int instance = 1; instance < p_base->number_particles_sector[1]; instance+=2) 
		quantum_number(0, p_base->number_particles_sector[0] + instance) = boundary_0 + instance;
	for (int instance = 0; instance < p_base->number_particles_sector[1]; instance+=2) 
		quantum_number(0, p_base->number_particles_sector[0] + instance) = - boundary_0 - (instance+1);
	// fill other sectors
	for (int string_type=1; string_type < p_base->numberTypes(); ++string_type) {
		for (int instance=1-p_base->number_slots_sector[string_type+1]%2; instance < p_base->number_particles_sector[string_type+1]; instance+=2) 
			quantum_number(string_type, instance) =  instance;	
		for (int instance=p_base->number_slots_sector[string_type+1]%2; instance < p_base->number_particles_sector[string_type+1]; instance+=2) 
			quantum_number(string_type, instance) =  - instance-1;
	}
	// sort the quantum numbers (they should be about the right order already)
	// the modifications below should not affect this sort order.
	// sorting order: 
	// 0, -2, 2, -4, 4, -6, 6, etc.
	// -1, 1, -3, 3, -5, 5, etc.
	sort(quantum_number, quantumNumberLessThan);
}



/** given ID constructor (default: base ground state (id 0); use id -1 (NO_ID) to *not* set id) **/
::State::State (const Base* const the_base, const long long int state_id)
	: 	p_chain(the_base->p_chain), p_base(the_base), quantum_number(the_base), 
		rapidity(the_base), convergence(NO_CONVERGENCE),
		its_lnnorm(NOT_CALCULATED), its_energy(NOT_CALCULATED), iterations(0), newton_iterations(0), 
		its_shifts(0), its_id(state_id), its_mode(NO_MODE), its_symmetry(NOT_SET)
{
	const char* here = "::State::State";
	if (state_id >= 0) setId(state_id);
	else if (NO_ID==state_id) its_id = NO_ID;
	else throw Exception (here, exc_InvalidID);
	//setFreeRapidities();
}

/** given shifts **/
::State::State (const Base* const the_base, const vector< Young >& shift)
	: 	p_chain(the_base->p_chain), p_base(the_base), quantum_number(the_base), 
		rapidity(the_base), convergence(NO_CONVERGENCE),
		its_lnnorm(NOT_CALCULATED), its_energy(NOT_CALCULATED), iterations(0), newton_iterations(0),
		its_shifts(), its_id(NO_ID), its_mode(NO_MODE), its_symmetry(NOT_SET)
{
	setShifts(shift);
}


/** given quantum numbers **/
::State::State (const Base* const on_base, const Strip<int>& quantum)
	: 	p_chain(on_base->p_chain), p_base(on_base), quantum_number(quantum), 
		rapidity(on_base), convergence(NO_CONVERGENCE),
		its_lnnorm(NOT_CALCULATED), its_energy(NOT_CALCULATED), iterations(0), newton_iterations(0), 
		its_shifts(), its_id(NO_ID), its_mode(NO_MODE), its_symmetry(NOT_SET)
{ 	
	// quantum numbers should be sorted for id() to work.
	sort(quantum_number, quantumNumberLessThan); 
}



/** copy constructor **/
::State::State (const State& original)
	: 	p_chain(original.p_chain), p_base(original.p_base), quantum_number(original.quantum_number), 
		rapidity(original.rapidity), convergence(original.convergence) ,
		its_lnnorm(original.its_lnnorm), its_energy(original.its_energy), 
		iterations(original.iterations), newton_iterations(original.newton_iterations),
		its_shifts(original.its_shifts), its_id(original.its_id), its_mode(original.its_mode), its_symmetry(original.its_symmetry)
{ }





/** assignment **/
/*
::State& ::State::operator= (const State& rhs)
{
	if (this != &rhs) {
		chain = rhs.chain;
		base = rhs.base; //this is problematic with references...
		quantum_number = rhs.quantum_number;
		rapidity = rhs.rapidity;
		convergence = rhs.convergence;
		its_lnnorm = rhs.its_lnnorm;
		its_energy= rhs.its_energy;
		iterations = rhs.iterations;
		newton_iterations = rhs.newton_iterations;
		its_shifts = rhs.its_shifts;
	}
	return *this;
}
*/


/** set quantum numbers from shifts **/
void ::State::setShift (const int sector, const Young& the_shift)
{
	// drat, we'll have to copy anyway.
	Young shift = the_shift;
	
	// first, calculate state_id==0 quantum numbers	in the right order to be able to apply the shifts 
	Strip<int> ground_quantum_number (p_base); 
	setGroundQuantumNumbers (p_base, ground_quantum_number);
	// keep track of easy id changes , if id and shifts both set
	if (0==sector && NO_ID != its_id && its_shifts.size() && (shift.height == its_shifts[0].height) && shift.height) {
		bool equal = true;
		for (int i=1; i < shift.height; ++i ) 
			if (shift.readAt(i) != its_shifts[0].readAt(i)) {
				equal = false;
				break;
			}
		if (equal) its_id += shift.readAt(0) - its_shifts[0].readAt(0);
		else its_id = NO_ID;
	}
	else its_id = NO_ID;	

	// if shifts field was set, keep it up to date.
	if (its_shifts.size()) its_shifts[sector] = shift;
	// transpose for sector 0
	if (0==sector) shift.transpose();

	// string type for this sector (sectors[0] and [1] are both type 0)
	int type = (sector<2)?0:sector-1;
	// apply shifts to quantum numbers, starting from the outside
	for (int particle=0; particle < p_base->number_particles_sector[sector]; ++particle) {
		int index = p_base->number_particles_sector[sector]-1-particle;
		// sector[1] starts after sector[0]
		if (1==sector) index += p_base->number_particles_sector[0]; 
		int shift_distance = (shift.readAt(particle)+(ground_quantum_number(type, index)>=0) )/2;
		int shift_sign = ((shift.readAt(particle)%2)?-1:1);
		quantum_number(type, index) = shift_sign * ground_quantum_number(type, index) + ((ground_quantum_number(type, index)<0)?-2:2)*shift_sign*shift_distance;  
	}

	// indicate uncalculated fields
	convergence = NO_CONVERGENCE;	
	its_lnnorm = NOT_CALCULATED;
	its_energy = NOT_CALCULATED;	
	iterations = 0;
	newton_iterations = 0;
	its_mode = NO_MODE;
	its_symmetry= NOT_SET;
}


/** set quantum numbers from shifts, return square diff with previous quantum numbers **/
// TODO: code copy-pasted from above...
int ::State::setShifts (const vector<Young>& the_shift)
{
	// drat, we'll have to copy anyway. there must be a way to avoid this.
	vector<Young> shift = the_shift;
	
	const char* here = "::State::setShifts";
	if (shift.size() != p_base->numberSectors()) throw Exception (here, exc_BadShift);

	/// check for acceptable Young tableaux (TODO: should be in Young)
	for (int i=0; i<shift.size(); ++i) {
		int last_width = shift[i].width;
		for (int j=0; j < shift[i].height; ++j) 
			if (shift[i].readAt(j) > last_width) throw Exception (here, exc_BadShift);
			else last_width = shift[i].readAt(j);
	}
	///
	
	its_shifts = shift;
	// first sector given as hole, not particle, shifts. 
	shift[0].transpose();
	// keep the old quantum numbers to calculate the difference at the end
	Strip<int> old_quantum_number (quantum_number);
	// first, calculate state_id==0 quantum numbers	in the right order to be able to apply the shifts 
	setGroundQuantumNumbers(p_base, quantum_number);
	// then, shift the quantum numbers according to shifts 
	for (int i=0; i < shift.size(); ++i) {
		// no shift: do nothing
		if (!shift[i].height || !shift[i].rows.size()) continue;
		// string type for this sector (sectors[0] and [1] are both type 0)
		int type = (i<2)?0:i-1;
		// apply shifts to quantum numbers, starting from the outside
		for (int particle=0; particle < p_base->number_particles_sector[i]; ++particle) {
			int index = p_base->number_particles_sector[i]-1-particle;
			// sector[1] starts after sector[0]
			if (1==i) index += p_base->number_particles_sector[0]; 
			int shift_distance = (shift[i].readAt(particle)+(quantum_number(type, index)>=0) )/2;
			// can't use sgn() because it returns 0 for 0 
			int shift_sign = ((shift[i].readAt(particle)%2)?-1:1);
			quantum_number(type, index) = shift_sign * quantum_number(type, index) + ((quantum_number(type, index)<0)?-2:2)*shift_sign*shift_distance;  
		}
	}
	// indicate uncalculated fields
	convergence = NO_CONVERGENCE;	
	its_lnnorm = NOT_CALCULATED;
	its_energy = NOT_CALCULATED;	
	iterations = 0;
	newton_iterations = 0;
	its_id = NO_ID;	
	its_mode = NO_MODE;
	its_symmetry = NOT_SET;
	// calculate distance from last state and return
	int dist = 0;
	for (int i=0; i< quantum_number.numberElements(); ++i) 
		dist += sq((quantum_number.element(i) - old_quantum_number.element(i))/2);
	return dist;
}



vector<Young> State::calculateShifts (void) const
{
	int number_sectors = p_base->numberSectors();

	// first, calculate state_id==0 quantum numbers	in the right order to be able to apply the shifts 
	Strip<int> ground_quantum_number (p_base);		
	setGroundQuantumNumbers(p_base, ground_quantum_number);
	
	// calculate young ids
	// the quantum numbers are assumed to be sorted!

	// sector[0]
	its_shifts.clear();
	vector<int> shift_0(p_base->number_particles_sector[0]);
	for (int instance=0; instance < p_base->number_particles_sector[0]; ++instance) {
		int index=p_base->number_particles_sector[0]-instance-1;
		shift_0[instance] = abs(quantum_number(0, index)) - abs(ground_quantum_number(0, index));
		if ((quantum_number(0, index)<0) != (ground_quantum_number(0, index)<0))
			shift_0[instance] += (quantum_number(0, index)<0)?-1:1;
	}
	its_shifts.push_back(Young(p_base->numberHoles(), p_base->number_particles_sector[0], shift_0));
	// other sectors
	for (int string_type=0; string_type < number_sectors-1; ++string_type) {
		if (!p_base->number_particles_sector[string_type+1]) continue;
		vector<int> shift(p_base->number_particles_sector[string_type+1]);
		
		for (int instance=0; instance < p_base->number_particles_sector[string_type+1]; ++instance) {
			int index=p_base->number_particles_sector[string_type+1]-instance-1;
			if (0==string_type) index+= p_base->number_particles_sector[0];
		///
			// 'difference' with ground quantum numbers
			// where distance is defined according to the sequence
			// 0, -2, 2, -4, 4, -6, 6, etc.
			// -1, 1, -3, 3, -5, 5, etc.
			shift[instance] = abs(quantum_number(string_type, index)) - abs(ground_quantum_number(string_type, index));
			// if the signs differ, we are one off in distance
			if ((quantum_number(string_type, index)<0) != (ground_quantum_number(string_type, index)<0))
				shift[instance] += (quantum_number(string_type, index)<0)?-1:1;
		///
		}
		its_shifts.push_back(Young(p_base->number_slots_sector[string_type+1]-p_base->number_particles_sector[string_type+1], p_base->number_particles_sector[string_type+1], shift));
	}
	
	// first sector given as hole, not particle, shifts. 
	its_shifts[0].transpose();
	return its_shifts;
}

/** set quantum numbers **/
void ::State::setQuantumNumbers (const Strip<int> the_quantum_numbers)
{
	const char* here = "::State::setQuantumNumbers";
	
	// copy the quantum numbers
	quantum_number = the_quantum_numbers;
	
	// indicate uncalculated fields
	its_shifts.clear();
	convergence = NO_CONVERGENCE;	
	its_lnnorm = NOT_CALCULATED;
	its_energy = NOT_CALCULATED;	
	iterations = 0;
	newton_iterations = 0;
	its_id = NO_ID;	
	its_mode = NO_MODE;
	its_symmetry = NOT_SET;
}
	
	
/** check admissibility of state **/
bool State::admissible(void) const
{
	// check quantum numbers. a few states need to be excluded 
	// because the string hypothesis leads to non-finite results
	
	// * the quantum numbers are symmetrically configured and there exist even-string-length bound states.
	//   - there will be an even string with real rapidity zero;
	//   - this implies (on the LHS of the original bethe equations) a lambda_j that is +- I zeta/2,
	//     yielding a zero numerator or denominator on the left, where no infinities should occur.
	//     (only (1+something)^N -type divergencies, being cancelled out by lambda_j - lambda_k = i zeta + delta on the right.)
	//  Solutions in this class are actually proper Bethe states (given the correct limiting procedure, see notes). 
	//  However, their contribution to correlations is necessarily zero (see notes).
	
	// * we have symmetrically distributed quantum numbers and there is more than one type of odd string of the same parity 
	//   e.g. a three-string with positive parity as well as the regular one-plusses.
	//   - we can have more than one zero rapidity (really, only at zero). this is forbidden by the exclusion principle.
	// 	NOTE: actually, this latter class is solvable by deviation. Also, its solutions do contribute to correlations.
	// 	So it must only be excluded if we don't deviate. Therefore, these should be counted as 'deviated', not 'inadmissible'
			
	// check for a zero in an even-length state. assume ordered quantum numbers.
	for (int j=0; j < p_base->numberTypes(); ++j) {
		if (p_base->numberStringsOfType(j) && !quantum_number(j, 0) && !(p_chain->stringLength(j)%2)) 
			return !symmetric();
	}	
	return true;
}
	


/** get number of holes **/
int ::State::calculateNumberHoles (void) const
{
	const char* here = "State::calculateNumberHoles";
	
	int number_holes = p_base->numberRoots();
	for (int i = 0; i < p_base->numberStringsOfType(0); ++i) {
		if (abs(quantum_number(0, i)) < p_base->numberRoots()) --number_holes;
		// left boundary of ground sector must be counted when parities don't match
		if (  ( (p_base->numberRoots()%2) != (p_base->numberStringsOfType(0)%2) ) &&  ( quantum_number(0, i) == -p_base->numberRoots() )  ) --number_holes; 	
	}
		
	if (number_holes<0) throw Exception(here, exc_NegativeHoles); 
	return number_holes;
}





/** one iteration of Newton's method **/
void ::State::newton (REAL& newt_convergence)
{
	const REAL tfh_threshold = 10.0;	// this is probably not generic...
	
	// Newton's method: Numerical Recipes in C,  p381
	int sign_permutation = 1;
	vector<int> indx (rapidity.numberElements()); // used by backsub, given by decompose 
	vector<REAL> delta (rapidity.numberElements());
	newt_convergence = 0.0;
	
	// Gaudin's matrix is Jacobi's matrix for Bethe's equations 
	Square<REAL> jacobian = matrixGaudin(); 

	// delta is the right side of   -F in mx_J - vc_deltax == - vec_F
	for (int j=0; j < p_base->numberTypes(); ++j) 
	for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha) 
		delta[rapidity.index(j,alpha)] = betheZero(j, alpha);

	decomposeLU (jacobian, indx, sign_permutation);
	backsubLU (jacobian, indx, delta);
	
	for (int i=0; i < rapidity.numberElements(); ++i) {
		if (abs(delta[i]) > tfh_threshold) continue; // don't generate tanhs that are 1
		if (isNan(delta[i])) continue;
		setRapidityAtIndex (i, rapidity.element(i) + delta[i]); 
		newt_convergence += sq(delta[i]);
	}
	
	++newton_iterations;
	iterate();			// iterate once to check convergence
}
	
	
	
/** simply iterate until convergence is reached **/
bool ::State::solve (const int max_iter, const REAL precision) 
{
	while ( (convergence > precision) && (++iterations < max_iter) ) iterate(); 
	return (convergence <= precision);
}



/** iterate a few times, then extrapolate  **/
bool ::State::solveExtrapolate (const int max_iter, const REAL precision) 
{
	const REAL rapidity_threshold = 18.0; // TODO NON_GENERIC!
	const int number_extrapolate = 4;  
	REAL dummy_convergence=0.0;
	REAL last_convergence=NO_CONVERGENCE;
	Strip<REAL>* rapidity_store[number_extrapolate];
	for (int i=0; i<number_extrapolate; ++i) {
		rapidity_store[i] = new Strip<REAL>(p_base);
	}
	vector<REAL> xvalues (number_extrapolate);
	vector<REAL> yvalues (number_extrapolate);
	while ( (convergence > precision) && (iterations < max_iter) ) {
		for (int i=0; i<number_extrapolate; ++i){
			iterate();
			if (convergence <= precision) break;
			*(rapidity_store[i]) = rapidity; 		//copy rapidity strip
		}
		if (convergence <= precision) break;
		last_convergence=convergence;
		// extrapolate along a trial-and-error-deduced magic curve
		vector<REAL> yvalues (number_extrapolate);
		for (int i=0; i< number_extrapolate; ++i) xvalues[i] = 2.0/(1.0*(i+4));	// the magic curve
		for (int j=0; j < rapidity.numberElements(); ++j) {
			for (int i=0; i<number_extrapolate; ++i) yvalues[i] = (*rapidity_store[i]).element(j);
			polint (xvalues, yvalues, 0.0, rapidity.element(j), dummy_convergence);

			// don't extrapolate to unreasonable values, fall back to last iteration.
			if (!finite(rapidity.element(j)) || (abs(rapidity.element(j)) > rapidity_threshold)) 
				setRapidityAtIndex(j, (*rapidity_store[number_extrapolate-1]).element(j) );
		}
		// check convergence
		iterate (); 
		// if extrapolation didn't help..
		if (convergence > last_convergence) 
			for (int j=0; j < rapidity.numberElements(); ++j) setRapidityAtIndex (j, (*rapidity_store[number_extrapolate-1]).element(j) );
	}
	for (int i=0; i<number_extrapolate; ++i) {
		 delete rapidity_store[i];
	}
	return (convergence <= precision);
}


/** iterate, extrapolate, newton's method: the whole shebang **/
bool ::State::solve (void)
{
	REAL newt_convergence = NO_CONVERGENCE;
	int last_iter_newton = 0;
	// first solve up to some precision
	solveExtrapolate (max_iterations, newton_threshold);
	while ( (convergence > precision) && (iterations < max_iterations)  ) {
		REAL last_convergence = convergence;
		// try newton's method within a given fluctuation band 
		// for at most a given number of times
		REAL last_newt_convergence = 2.0*newt_convergence;  
		last_iter_newton = newton_iterations;
		while ((newt_convergence < newton_bandwidth * last_newt_convergence)
				&& newton_iterations < max_newton
				&& newton_iterations < last_iter_newton + newton_consecutive
				) {
			last_newt_convergence = newt_convergence;
			try { 
				newton (newt_convergence); 
			}
			catch (const char* msg) {
				// must be decomposeLU: singular matrix.
				// newton can't be used, see if iterating once helps
				its_lnnorm = NOT_CALCULATED;
				iterate ();
			}
			
			if (convergence<=precision) break; 
			else its_lnnorm = NOT_CALCULATED;
		}
		if (convergence<=precision) break;
		
		// if we're here, newton failed for now.
		its_lnnorm = NOT_CALCULATED;
		
		if (newton_iterations >= max_newton) {
			// we have reached the max. no other choice but to try interpolation
			// if we reach max_iter, we haven't converged.
			solveExtrapolate (max_iterations, precision);
			break;
		}
		else {
			// we're here, so we haven't reached the max.
			// however, newton is not an option right now.

			// go a given factor from the current convergence.
			// do only a limited number of iterations, then newton again.
			newton_threshold = max(newton_factor * convergence, precision);
			solveExtrapolate (iterations + extrapolate_consecutive, newton_threshold);
		}
	}
	return (convergence <= precision);
}

	

/** momentum (mode number), between 0 (inclusive) and number_sites (exclusive) **/
int ::State::calculateMode (void) const
{	
	int number_pos_parity=0;
	for (int string_type=0; string_type < p_chain->numberTypes(); ++string_type)
		if (p_chain->positiveParity(string_type))  number_pos_parity += p_base->numberStringsOfType(string_type);
	int sum_quantum_numbers=0;
	for (int q=0; q< numberRapidities(); ++q) sum_quantum_numbers += quantum_number.element(q);
	// quantum number == 2 I
	its_mode = (p_chain->length()*(number_pos_parity) + sum_quantum_numbers)/2;
	// bring to default domain (-pi, pi]
	while (its_mode >= p_chain->length()) its_mode -= p_chain->length();
	while (its_mode < 0) its_mode += p_chain->length();
	return its_mode;
}

bool State::calculateSymmetry(void) const
{
	// let's see if the distribution of quantum numbers is symmetric.
	/*
	// this should work for unordered quantum numbers
	for (int j=0; j < base.numberTypes(); ++j) { 
		int sum_quantum_numbers = 0, sum_cube_quantum_numbers = 0;
		for (int alpha=0; alpha < base.numberStringsOfType(j); ++alpha) {
			sum_quantum_numbers += quantum_number(j, alpha);
			sum_cube_quantum_numbers += quantum_number(j, alpha)*quantum_number(j, alpha)*quantum_number(j, alpha);
		}
		if (sum_quantum_numbers || sum_cube_quantum_numbers) return false;
	}
	return true;
	*/
	
	// faster: assume ordering 0 -1 1 -2 2 etc which we enforce anyway
	for (int j=0; j < p_base->numberTypes(); ++j) { 
		int number_j = p_base->numberStringsOfType(j);
		// central string not zero: not symmetric
		if (number_j%2 && quantum_number(j,0)) return false;
		for (int alpha=number_j%2; alpha < number_j; alpha+=2) {
			// if not symmetric, bail out immediately
			if (quantum_number(j, alpha+1) != -quantum_number(j, alpha)) return false;
		}
	}
	return true;
}		

