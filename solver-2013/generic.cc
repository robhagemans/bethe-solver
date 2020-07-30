#include "generic.h"




// monikers
const string mon_File = "v20105";
const string mon_Base = "base";
const string mon_Delta = "D";
const string mon_Length = "N";
const string mon_LeftDown = "M";
const string mon_RightDown = "R";
const string mon_Spinon = "s";

const string sep_Name = "_"; 
const string sep_Final = ".";
const string sep_Vector = "-";

const char* exc_Ferro = "ferromagnetic and XX0 chains not supported";		// generic contructors
const char* exc_InfiniteRaps = "only infinite rapidities for the isotropic chain";
const char* exc_NegativeAnisotropy = "negative or zero anisotropy"; 			//name()
const char* exc_BadName = "cannot parse the moniker";
const char* exc_ThresholdOne = "contribution threshold must be less than one";
const char* exc_NullPointer = "null pointer passed";
const char* exc_NotConverged = "bethe equations not converged";	//solve()
const char* exc_Deviated = "bethe strings deviated";	//solve()



/** generic construction **/
/** new chain **/

Chain* newChain (const double delta, const int chain_length, const int cutoff_types) 
{
	const string here = "newChain";
	if (delta > 1.0) return new Gap_Chain (delta, chain_length, cutoff_types);
	if (delta == 1.0) return new XXX_Chain (chain_length, cutoff_types);
	if ((delta < 1.0) && (delta > 0.0)) return new XXZ_Chain (delta, chain_length);
	throw Exception (here, exc_Ferro);	
}

/** new base **/

Base* newBase (const Chain* const p_chain, const BaseData& base_data) {
	const string here = "newBase";
	if (p_chain->delta() == 1.0) return new XXX_Base (p_chain, base_data);
// 	if (number_infinite_rapidities) throw Exception (here, exc_InfiniteRaps);	
	if (p_chain->delta() > 1.0) return new Gap_Base (p_chain, base_data);
	if ((p_chain->delta() < 1.0) && (p_chain->delta() > 0.0)) return new XXZ_Base (p_chain, base_data);
	throw Exception (here, exc_Ferro);	
}


Base* newBase (const Chain* const p_chain, const vector<int>& base_vec, const int number_spinons, const int number_infinite_rapidities)
{
	const string here = "newBase";
	if (p_chain->delta() == 1.0) return new XXX_Base (p_chain, base_vec, number_spinons, number_infinite_rapidities);
	if (number_infinite_rapidities) throw Exception (here, exc_InfiniteRaps);	
	if (p_chain->delta() > 1.0) return new Gap_Base (p_chain, base_vec, number_spinons);
	if ((p_chain->delta() < 1.0) && (p_chain->delta() > 0.0)) return new XXZ_Base (p_chain, base_vec, number_spinons);
	throw Exception (here, exc_Ferro);	
}

Base* newBase (const Chain* const p_chain, const int number_down, const vector<int>& base_vec, const int number_spinons, const int number_infinite_rapidities)
{
	const string here = "newBase";
	if (p_chain->delta() == 1.0) return new XXX_Base (p_chain, number_down, base_vec, number_spinons, number_infinite_rapidities);
	if (number_infinite_rapidities) throw Exception (here, exc_InfiniteRaps);	
	if (p_chain->delta() > 1.0) return new Gap_Base (p_chain, number_down, base_vec, number_spinons);
	if ((p_chain->delta() < 1.0) && (p_chain->delta() > 0.0)) return new XXZ_Base (p_chain, number_down, base_vec, number_spinons);
	throw Exception (here, exc_Ferro);	
}

Base* newGroundBase (const Chain* const p_chain, const int number_down)
{
	const string here = "newGroundBase";
	if (p_chain->delta() > 1.0) return new Gap_Base (p_chain, number_down);
	if (p_chain->delta() == 1.0) return new XXX_Base (p_chain, number_down);
	if ((p_chain->delta() < 1.0) && (p_chain->delta() > 0.0)) return new XXZ_Base (p_chain, number_down);
	throw Exception (here, exc_Ferro);	
}


/** new state **/

State* newGroundState (const Base* const p_base) 
{
	const string here = "newGroundState";
	if (p_base->p_chain->delta() > 1.0) return new Gap_State (p_base);
	if (p_base->p_chain->delta() == 1.0) return new XXX_State (p_base);
	if ((p_base->p_chain->delta() < 1.0) && (p_base->p_chain->delta() > 0.0)) return new XXZ_State (p_base);
	throw Exception (here, exc_Ferro);	
}

State* newState (const Base* const p_base, const long long int id) 
{
	const string here = "newState";
	if (p_base->p_chain->delta() > 1.0) return new Gap_State (p_base, id);
	if (p_base->p_chain->delta() == 1.0) return new XXX_State (p_base, id);
	if ((p_base->p_chain->delta() < 1.0) && (p_base->p_chain->delta() > 0.0)) return new XXZ_State (p_base, id);
	throw Exception (here, exc_Ferro);	
}

State* newState (const Base* const p_base, const vector<Young>& shifts) 
{
	const string here = "newState";
	if (p_base->p_chain->delta() > 1.0) return new Gap_State (p_base, shifts);
	if (p_base->p_chain->delta() == 1.0) return new XXX_State (p_base, shifts);
	if ((p_base->p_chain->delta() < 1.0) && (p_base->p_chain->delta() > 0.0)) return new XXZ_State (p_base, shifts);
	throw Exception (here, exc_Ferro);	
}


/** clone **/

// TODO: replace by State::clone()
State* copy (State* const p_state) 
{
	const string here = "State* copy";
	if (!p_state) throw Exception (here, exc_NullPointer);
	if (p_state->p_chain->delta() > 1.0) return new Gap_State ( *( (Gap_State*) p_state ));
	if (p_state->p_chain->delta() == 1.0) 	return new XXX_State ( *( (XXX_State*) p_state ));
	if ((p_state->p_chain->delta() < 1.0) && (p_state->p_chain->delta() > 0.0)) return new XXZ_State ( *( (XXZ_State*) p_state ));
	throw Exception (here, exc_Ferro);	
}



/** file naming conventions **/

/** parse names **/

string nameFromValue (const string moniker, const REAL value)
{
	stringstream namestream;
	namestream << sep_Name << moniker << value;
	return namestream.str();
}

string nameFromValue (const string moniker, const int value, const int width)
{
	stringstream namestream;
	namestream << sep_Name << moniker;
	if (width) {
		namestream.width (width);
		namestream.fill('0');
	}
	namestream << value;
	return namestream.str();
}

string nameFromValue (const string moniker, const vector<int>& value)
{
	stringstream namestream;
	namestream << sep_Name << moniker;
	for (int i=0; i<value.size(); ++i) namestream << sep_Vector << value[i];
	return namestream.str();
}


int intValueFromName(const string name, const string moniker)
{
	const string here = "intValueFromName";
	int length= moniker.length() + sep_Name.length();
	// find anisotropy delta
	int start = name.find(sep_Name + moniker, 0);
	int finish = name.find(sep_Name, start+length); // number ends at an underscore
	if (finish == string::npos) finish = name.find(sep_Final, start+length); // or at a dot
	if ((start==string::npos) || (finish==string::npos)) throw Exception(here, exc_BadName, moniker);
	return atoi(name.substr(start+length, finish-start-length).c_str());
}

REAL realValueFromName(const string name, const string moniker)
{
	const string here = "realValueFromName";
	int length = moniker.length() + sep_Name.length();
	// find anisotropy delta
	int start = name.find(sep_Name + moniker, 0);
	int finish = name.find(sep_Name, start+length); // number ends at an underscore
	if ((start==string::npos) || (finish==string::npos)) throw Exception(here, exc_BadName, moniker);
	return atof(name.substr(start+length, finish-start-length).c_str());
}

string valueFromName(const string name, const string moniker)
{
	const string here = "valueFromName";
	int length = moniker.length() + sep_Name.length();
	// find anisotropy delta
	int start = name.find(sep_Name + moniker, 0);
	int finish = name.find(sep_Name, start+length); // number ends at an underscore
	if (finish == string::npos) finish = name.find(sep_Final, start+length); // or at a dot
	if (finish == string::npos) finish = name.length(); // or at the end
	if ((start==string::npos) || (finish==string::npos)) throw Exception(here, exc_BadName, moniker);
	return name.substr(start+length, finish-start-length);
}





/** base name **/
string name(const BaseData& base_data) 
{
	stringstream the_name;
	//the_name.precision(3);
	the_name << nameFromValue(mon_RightDown, base_data.number_magnons)  << nameFromValue(mon_Spinon, base_data.number_spinons)	<< nameFromValue(mon_Base, base_data.its_structure);
	return the_name.str();
}

string name(const Base* const p_base) 
{
	stringstream the_name;
	//the_name.precision(3);
	the_name << nameFromValue(mon_RightDown, p_base->numberDown())  << nameFromValue(mon_Spinon, p_base->numberSpinons())	<< nameFromValue(mon_Base, p_base->structure());
	return the_name.str();
}

/** chain name **/
string name(const Chain* const p_chain, const int number_down)
{
	const char* here = "name";
	stringstream description;
	description << sep_Name << mon_File;
		
	if (p_chain->delta() <= 0) throw Exception (here, exc_NegativeAnisotropy);
	if (p_chain->delta() == 1.0) description << "_xxx";
	if (p_chain->delta() > 1.0) description << "_xxz-gap";
	if (p_chain->delta() < 1.0) description << "_xxz";
	
	description << nameFromValue(mon_Delta, p_chain->delta());
	description << nameFromValue(mon_Length, p_chain->length());
	if (number_down) description << nameFromValue (mon_LeftDown, number_down, fieldWidth(p_chain->length()) ); // all fields from 1 to N/2 equally wide
	return description.str();
}

// convert a name into a chain
Chain* newChain (const string name, const int cutoff_types) 
{	
	REAL delta = realValueFromName (name, mon_Delta);
	int chain_length = intValueFromName (name, mon_Length);
	return newChain (delta, chain_length, cutoff_types); 
}


// convert a name into a base
Base* newBase (const string name, Chain*& p_chain)
{
	const string here = "newBase";
	// if no chain is defined, get our own (note: passing a literal zero for Chain* guarantees trouble, you need to keep the result!)
	
	int right_number_down = intValueFromName(name, mon_RightDown); // number of down spins (right side) incl. infinite
	int number_spinons = intValueFromName(name, mon_Spinon);
	vector<string> base_elements = explode(valueFromName(name, mon_Base), sep_Vector);
	// if necessary, create the chain
	if (!p_chain) p_chain = newChain(name, base_elements.size()+1);
	// find the number of infinite rapidities & create a base vector
	int number_down_found=0;
	vector<int> base_vec (base_elements.size());
	for (int i=0; i<base_elements.size(); ++i) {
		base_vec[i] = atoi(base_elements[i].c_str());
		// this is only necessary due to the silly way I implemented the constructors.
		number_down_found += base_vec[i] * p_chain->stringLength(i);
	}
	/// TODO: implement a contructor that takes base_vec on its word as a boy scout and calculates the number of infinite raps by itself from number_down.
	return newBase(p_chain, right_number_down, base_vec, number_spinons, right_number_down - number_down_found);
}



/// BaseData constructor?
BaseData readBaseData (const string name) {
	int right_number_down = intValueFromName(name, mon_RightDown); // number of down spins (right side) incl. infinite
	int number_spinons = intValueFromName(name, mon_Spinon);
	vector<string> base_elements = explode(valueFromName(name, mon_Base), sep_Vector);
	int number_down_found=0;
	vector<int> base_vec (base_elements.size());
	for (int i=0; i<base_elements.size(); ++i) {
		base_vec[i] = atoi(base_elements[i].c_str());
		// this is only necessary due to the silly way I implemented the constructors.
// 		number_down_found += base_vec[i] * p_chain->stringLength(i);
	}
	return BaseData(right_number_down, base_vec, NOT_SET, number_spinons);
}



/** scan through bases, particle/spinon ordered **/
vector<BaseData> allNewBases (
	const Chain* const p_chain, const int number_down, 
	const int max_string_length,
	const int uptoinc_number_particles, const int uptoinc_number_spinons,
	const int max_infinite
) 
{
	vector<BaseData> all_bases;
	vector<int> basevec (1, 0); // set dummy first element
	int number_particles = 0;
	int last_type_increased = 0;
	bool limits_exceeded = false;
	//int max_infinite = (p_chain->delta() == 1.0) ?2 :0;
	while (number_particles <= min(uptoinc_number_particles, number_down)) {
		limits_exceeded = false;
		// every higher excitation contributes itself as a particle 
		int number_higher_particles = 0;
		for (int i=1; i < basevec.size(); ++i) number_higher_particles += basevec[i] ;
		if (number_higher_particles > number_particles) limits_exceeded=true; // all particles must be higher particles
		// we must fit the available number of particles exactly, so we need an even number of spinons+holes remaining.
		else if (number_higher_particles == number_particles) { 
			for (int number_spinons = 0; number_spinons <= min(uptoinc_number_spinons, number_down); ++number_spinons ) {
				for (int number_infinite=0; number_infinite <= max_infinite; ++number_infinite){
					try {
						// try to make a base, will throw if it can't
						// NOTE: this constructor calculates the number of type-1 particles from the other info
						Base* p_base = newBase(p_chain, number_down, basevec, number_spinons, number_infinite);
						// but keep only the base data (as calculated by the constructor!)
						BaseData base_data = *p_base;
						all_bases.push_back(base_data);
					}
					catch (Exception exc) {
						if (exc.error == exc_BlackburnLancashire
							|| exc.error == exc_TooMany) {
							// too many holes for this base , i.e. not enough space for spinons. ignore.
						}
						else throw;
					}

				}
			}
		}
		// find out where to increase 
		int number_excited;
		int type_to_increase = 1;
		if (limits_exceeded) type_to_increase = last_type_increased +1;
		if ((type_to_increase >= p_chain->numberTypes()) || (p_chain->stringLength(type_to_increase) > max_string_length )) {
			// we're out of bases for this number of particles. increase and reset p_base->
			++number_particles;
			basevec.clear();
			basevec.push_back(0); // set the dummy first element
			type_to_increase = 0; // increase only the dummy element, aka start with the trivial p_base->
		}
		// increase.
		for (int type=0; type < type_to_increase; ++type) basevec[type] = 0;
		if (type_to_increase >= basevec.size()) basevec.push_back(0);
		++basevec[type_to_increase];
		last_type_increased = type_to_increase;
	}	
	return all_bases;
}




/** generic solve **/
State* solve(
	State* const p_state_in,
	const REAL deviation_threshold,
	const bool deviate
)
{
	const char* here = "::solve()";
	// don't clone, just copy 
	// cloning only necessary when deviating
	// TODO: by separating id info from solving and roots, we could 
	// do this as fast and much more elegantly.
	State* p_state = p_state_in;
	
	// solve the Bethe Equations up to convergence
	if (p_state->solve()) { 
		// our state has converged
		REAL deviation = p_state->stringDeviation();
		// our core business: calculate form factor
		if (deviation < deviation_threshold) return p_state; 
		else if (deviate && p_state->p_chain->delta()==1.0) {
			XXXDeviatedState* p_dev_state = new XXXDeviatedState ( *p_state );
			if (!p_dev_state->solve()) {
				// we haven't converged! (counts as deviated)
				stringstream message;
				message << "deviance not converged, convergence "<<p_dev_state->convergence<< " threshold "<< State::precision;
				// clean up
				delete p_dev_state;
				throw Exception(here, exc_Deviated, message.str());
			}
			else return p_dev_state;
		}
		else {
			stringstream message;
			message << "non-deviable, deviation threshold "<<deviation_threshold;
			throw Exception(here, exc_Deviated, message.str());
		}
	}
	else {
		stringstream message;
		message << "convergence "<<p_state->convergence<< " threshold "<< State::precision;
		throw Exception(here, exc_NotConverged, message.str());
	}
}
