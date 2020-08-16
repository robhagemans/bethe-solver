
#include "base.h"
#include "state.h"
#include "generic.h"
#include "quantity.h"

const string Longitudinal::ff_name = "Szff";
const string Longitudinal::sf_name = "Szz";
const string Transverse::ff_name = "Smff";
const string Transverse::sf_name = "Smp";
const string SPlusMinus::ff_name = "Spff";
const string SPlusMinus::sf_name = "Spm";


// convert a name into a Quantity
Quantity* newQuantityFromFileName (const string name, State* const p_ground_state)
{
	const char* here = "newQuantity";
	string function_name = name.substr(0, name.find (sep_Name)); // function name is first token
	Quantity* p_quantity = newQuantityFromName (function_name, p_ground_state);
	if (!p_quantity) throw Exception(here, exc_BadName, function_name);
	return p_quantity;
}




/** magnetic fields for which the given number of down spins is the ground state **/

vector<REAL> magneticFields (const Chain* const p_chain)
{
	int iterations, newton_iterations;

	vector<REAL> energy (p_chain->length()/2+1);
	vector<REAL> field (p_chain->length()/2+1);
	REAL diff_energy=0.0, last_diff_energy=0.0;
	field[0] = 2.0;
	field[p_chain->length()/2] = 0.0;

	Base* p_base = newGroundBase (p_chain, 0);
	State* p_state = newGroundState (p_base);
	p_state->solve();
	energy[0] = p_state->energy();
	delete p_base, p_state;

	for (int number_down=0; number_down < p_chain->length()/2; ++number_down) {
		// construct a subspace ground base & ground state
		Base* p_base = newGroundBase (p_chain, number_down+1);
		State* p_state = newGroundState  (p_base);
		p_state->solve();
		energy[number_down+1] = p_state->energy();
		if (number_down) field[number_down] = (energy[number_down-1]-energy[number_down+1])*0.5;
		delete p_base;
		delete p_state;
	}

	return field;
}

REAL magneticField (const Chain* const p_chain, const int number_down)
{
	if (number_down==0) return 2.0;
	if (number_down==p_chain->length()/2) return 0.0;

	int iterations, newton_iterations;

	Base* p_base_hi = newGroundBase (p_chain, number_down+1);
	State* p_state_hi = newGroundState (p_base_hi);
	p_state_hi->solve();

	Base* p_base_lo = newGroundBase  (p_chain, number_down-1);
	State* p_state_lo = newGroundState  (p_base_lo);
	p_state_lo->solve();

	REAL field = (p_state_lo->energy()-p_state_hi->energy())*0.5;
	delete p_base_hi;
	delete p_state_hi;
	delete p_base_lo;
	delete p_state_lo;
	return field;
}

