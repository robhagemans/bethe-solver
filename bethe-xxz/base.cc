#include "base.h"


const char* exc_UpsideDown = "more down spins than half the chain length";
const char* exc_SadMouse = "too few holes for this base";
const char* exc_BlackburnLancashire = "too many holes for this base";
const char* exc_NumberTypes = "too many types";
const char* exc_SwissCheese = "more holes than down spins";
const char* exc_TooMany = "more particles than down spins";
const char* exc_Limits = "base exceeds limits: more particles than slots for some sector";
const char* exc_HolesNotSet = "incomplete construction: number_holes not set";
const char* exc_IDTooHigh = "ID too high";
const char* exc_Infinite = "no infinite rapidities in anisotropic chain";

// sentinel value for naturals
const int NOT_SET = -666;



/** class BaseData **/

/** default empty constructor **/
BaseData::BaseData (void): its_structure(), number_magnons(0), number_holes(0) {};

/** set structure **/
BaseData::BaseData (const int the_number_down, const vector<int>& the_structure, const int the_number_holes, const int the_number_spinons):
	its_structure(the_structure), number_holes(the_number_holes), number_magnons(the_number_down), number_spinons(the_number_spinons)
{
	crop();
}

/** copy **/
BaseData::BaseData (const BaseData& original):
	its_structure(original.its_structure), number_holes(original.number_holes), number_magnons(original.number_magnons) , number_spinons(original.number_spinons)
{
	crop();
}

/** assign **/
BaseData& BaseData::operator= (const BaseData& rhs)
{
	if (this != &rhs) {
		its_structure = rhs.its_structure;
		number_holes = rhs.number_holes;
		number_spinons = rhs.number_spinons;
		number_magnons = rhs.number_magnons;
	}
	return *this;
}

/** remove trailing zeroes, add zero to empty vector **/
void BaseData::crop(void)
{
	// remove trailing zeroes from the structure vector: they spoil numberTypes.
	// leave at least one entry, for base 0 has number types 1.
	if (its_structure.size()) while (!its_structure.back() && its_structure.size()>1 ) its_structure.pop_back();
	else its_structure.push_back(0);
}

/** number of degrees of freedom **/
int BaseData::numberFreedoms (void) const
{
	// a freedom is one of:
	//   - a hole in the ground interval
	//   - a spinon
	//   - a string (1-, 2, etc) regardless of its length
	// notably, an infinite rapidity is NOT a freedom and a 2-string is only 1 freedom

/// FIXME: we do it totally different now
/// this is basically the number of holes. for XXX, but it's a tweak.
/// it should be implemented better, but I'm sleepy.
	int number_freedom = number_spinons; //+ number_holes;
	for (int i=1; i<its_structure. size(); ++i) number_freedom += its_structure[i]*(i+1);
	return number_freedom;
}




/** class Base **/

/** construct from baseData **/
Base::Base (const Chain* const on_chain, const BaseData& base_data) : BaseData(base_data), p_chain(on_chain), number_slots_sector(0), lim_id (0)
{
	const char* here = "Base::Base:";

	if (its_structure.size() > p_chain->numberTypes()) throw Exception(here, exc_NumberTypes);
	// set the length of the other vectors
	number_slots_sector.resize(its_structure.size()+1);
	lim_id.resize(its_structure.size()+2);
	// this constructor believes number down and calculates the infinite raps (if any) from it.
}



/** given-base constructor  **/
Base::Base (const Chain* const on_chain, const vector<int>& structure, const int new_number_spinons)
	: BaseData (0, structure, NOT_SET, new_number_spinons), p_chain(on_chain), number_slots_sector(0), lim_id (0)
{
	const char* here = "Base::Base:";

	if (its_structure.size() > p_chain->numberTypes()) throw Exception(here, exc_NumberTypes);
	// set the length of the other vectors
	number_slots_sector.resize(its_structure.size()+1);
	lim_id.resize(its_structure.size()+2);
	// calculate number down
	for (int i=0; i< numberTypes(); ++i) number_magnons += numberStringsOfType(i) * p_chain->stringLength(i);
	// setNumberHoles/Spinons is not called. this needs calculateNumberSlots(), which is pure virtual--> trouble.
}

/** given-base constructor, ignore first element of structure and calculate it from number down **/
Base::Base (const Chain* const on_chain, const int the_number_roots, const vector<int>& structure, const int new_number_spinons)
	: BaseData (the_number_roots, structure, NOT_SET, new_number_spinons), p_chain(on_chain), number_slots_sector(0), lim_id (0)
{
	const char* here = "Base::Base:";

	if (its_structure.size() > p_chain->numberTypes()) throw Exception(here, exc_NumberTypes);
	// empty vector gets modified to ground base. --- now in crop()
	// set the length of the other vectors
	number_slots_sector.resize(its_structure.size()+1);
	lim_id.resize(its_structure.size()+2);

	// calculate first element from given number of down spins
	its_structure[0] = the_number_roots;
	for (int i=1; i< its_structure.size(); ++i) {
		its_structure[0] -= its_structure[i] * p_chain->stringLength(i);
	}
	if (its_structure[0]<0) throw Exception(here, exc_TooMany);
}


/** ground-base constructor **/
Base::Base (const Chain* const on_chain, const int the_number_roots)
	:	BaseData(the_number_roots, vector<int> (1, the_number_roots), 0, 0), p_chain(on_chain), number_slots_sector(its_structure.size()+1), lim_id (its_structure.size() + 2)
{}


/** copy constructor **/
Base::Base (const Base& original)
	: 	BaseData(original),
		p_chain(original.p_chain),
		number_slots_sector(original.number_slots_sector),
		number_particles_sector(original.number_particles_sector),
		lim_id(original.lim_id)
{}


/** assign **/
Base& Base::operator= (const Base& rhs)
{
	if (this != &rhs) {
		(BaseData) *this = (BaseData) rhs;
		p_chain = rhs.p_chain;
		number_slots_sector = rhs.number_slots_sector;
		number_particles_sector = rhs.number_particles_sector;
		lim_id= rhs.lim_id;
	}
	return *this;
}


/** set structure **/
Base& Base::setStructure (const vector<int>& structure, const int new_number_spinons)
{
	const char* here = "Base::setStructure";
	its_structure.clear();
	its_structure = structure;
	crop();

	if (its_structure.size() > p_chain->numberTypes()) throw Exception(here, exc_NumberTypes);
	// empty vector gets modified to ground base.
	if (!its_structure.size()) its_structure.push_back(0);
	// set the length of the other vectors
	number_slots_sector.resize(its_structure.size()+1);
	lim_id.resize(its_structure.size()+2);

	// calculate number of down spins
	number_magnons=0;
	for (int i=0; i< its_structure.size(); ++i)
		number_magnons += its_structure[i] * p_chain->stringLength(i);

	setNumberSpinons(new_number_spinons);
	calculateNumberSlots();
	checkConsistency();
	return *this;
}


/** set number holes **/
Base& Base::setNumberHoles (const int new_number_holes)
{
	// a hole is an unfilled slot in the ground sector.
	// should the ground sector be larger than the number of particles of string type 1+,
	// as can happen in the case of parity mismatch between that number and the number of slots for string type 1+,
	// (depending on the constant PLUS_OR_MINUS_ONE_TO_SYMMETRISE) there is always a hole in the ground sector.

	const char* here = "Base::setNumberHoles";
	if (new_number_holes > numberDown() ) throw Exception(here, exc_SwissCheese);
	if (new_number_holes > maxNumberHoles()) throw Exception (here, exc_BlackburnLancashire);
	if (new_number_holes < minNumberHoles()) throw Exception (here, exc_SadMouse);
	number_holes = new_number_holes;
	number_spinons = number_holes - minNumberHoles();
	calculateNumberParticles ();
	calculateLimId ();
	return *this;
}


/** set number spinons **/
Base& Base::setNumberSpinons (const int new_number_spinons)
{
	const char* here = "Base::setNumberSpinons";
	number_spinons = new_number_spinons;
	number_holes = minNumberHoles() + number_spinons;

	if (number_holes > maxNumberHoles()) throw Exception (here, exc_BlackburnLancashire);
	if (number_holes < minNumberHoles()) throw Exception (here, exc_SadMouse);

	calculateNumberParticles ();
	calculateLimId ();
	return *this;
}



/** consistency checks **/
void Base::checkConsistency (void) {
	const char* here = "Base::checkConsistency";

	if (2*numberDown() > p_chain->length()) throw Exception (here, exc_UpsideDown);
	if (number_particles_sector[0] > number_slots_sector[0]) throw Exception (here, exc_SadMouse);
	for (int i=1; i< numberSectors(); ++i)
		if (number_particles_sector[i] > number_slots_sector[i]) throw Exception (here, exc_BlackburnLancashire);
}


/** calculate number of slots per sector **/
void Base::calculateNumberSlots (void)
{
	const char* here = "Base::calculateNumberSlots";

	number_slots_sector = limQuantumNumbers(); // copy

	// even number particles --> odd quantum numbers --> even number of slots
	// and vice versa.
	for (int string_type=0; string_type < numberTypes(); ++string_type)
		if ((numberStringsOfType(string_type)%2) != (number_slots_sector[string_type]%2)) --number_slots_sector[string_type];
	// consistency check
	for (int i=0; i < numberTypes(); ++i)
		if ( numberStringsOfType(i) > number_slots_sector[i] ) throw Exception (here, exc_Limits);

	// now interval_width = number slots per string type.

	// in the ground state, if number_roots is even, I s are half-integers
	// in other words, quantum_numbers of the ground state are odd.
	// in that case, the number of slots in the ground *sector* is the number of odd integers below 2 *max_quantum (i.e. below 4*I_max)
	// which equals lim_quantum if even, or lim_quantum-1 if odd.

	int ground_interval_width = min(number_magnons, numberStringsOfType(0)); 	// this choice seems to be compatible with J-S's code

	// keep the ground interval symmetric: make sure it has the same even-odd parity as the number of slots it fits into.
	if (ground_interval_width%2 != number_slots_sector[0]%2) ground_interval_width += PLUS_OR_MINUS_ONE_TO_SYMMETRISE;

	// first sector: number slots
	number_slots_sector.insert (number_slots_sector.begin(), ground_interval_width);
	number_slots_sector[1] -= number_slots_sector[0]; // sector [1] is what's left after removing sector [0].

}


/** calculate number of particles per sector **/
void Base::calculateNumberParticles (void)
{
	const char* here = "Base::calculateNumberParticles";
	// needs number_holes to be set
	if (number_holes == NOT_SET) throw Exception (here, exc_HolesNotSet);
	// copy. this is stupid.
	number_particles_sector = its_structure;
	// where there's a hole, there's no particle.
	number_particles_sector.insert (number_particles_sector.begin(), number_slots_sector[0] - number_holes);
	// what's in sector 0 is not in sector 1.
	number_particles_sector[1] -= number_particles_sector[0];
}



/** calculate the lowest integer that is not a valid id number **/
void Base::calculateLimId (void)
{
	const char* here = "Base::calculateLimId";

	lim_id[0] = 1;
	for (int i=1; i<= numberSectors(); ++i) {
		lim_id[i] = lim_id[i-1] * choose (number_slots_sector[i-1], number_particles_sector[i-1]);
		if (lim_id[i] < 0) throw Exception (here, exc_Overflow);
	}
	while (!lim_id.back()) lim_id.pop_back();
}





/** number of slots per string type **/
vector<int> Base::numberSlotsPerStringType (void) const
{
	vector<int> number_slots_type = number_slots_sector;
	number_slots_type[1] += number_slots_type[0];
	number_slots_type.erase(number_slots_type.begin());
	return number_slots_type;
}


/** get the id number of a set of shifts **/
long long int Base::id (const vector<Young>& shifts) const
{
	// get state id from young ids.
	long long int the_id = 0;
	for (int i = numberSectors()-1; i >=0; --i) {
		// number of young tableaux
		the_id *= choose (number_slots_sector[i], number_particles_sector[i]);
		the_id += shifts[i].id();
	}
	return the_id;
}

// get the shifts belonging to a certain id number
vector<Young> Base::shifts(const long long int id) const
{
	const char* here = "Base::shifts";

	// calculate shifts from state_id
	int number_sectors = numberSectors();
	if (id >= lim_id[number_sectors]) throw Exception (here, exc_IDTooHigh);

	// get a list of young_ids for each sector from the state id
	vector<long long int> young_id (number_sectors);
	long long int current_id = id;
	for (int i = number_sectors-1; i >= 0; --i) {
		young_id[i] = (long long int)(floor( 1.0 * current_id / lim_id[i] ));
		current_id -= young_id[i] * lim_id[i];
	}
	vector<Young> the_shifts;
	for (int i=0; i < number_sectors; ++i) {
		// calculate shifts: get rows in young tableaux from tableau id.
		Young particle_shift;
		if (0==i) {
			// for sector 0, move holes, then calculate particle shift.
			particle_shift.setId( number_particles_sector[0], number_holes, young_id[i] );
		}
		else particle_shift.setId( number_slots_sector[i] - number_particles_sector[i], number_particles_sector[i], young_id[i] );
		the_shifts.push_back(particle_shift);
	}
	return the_shifts;
}
