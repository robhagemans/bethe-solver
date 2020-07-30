#ifndef BASE_H
#define BASE_H

#include <vector>
#include "bethe.h"
#include "chain.h"
#include "young.h"

extern const char* exc_UpsideDown;
extern const char* exc_SadMouse;
extern const char* exc_BlackburnLancashire;
extern const char* exc_NumberTypes;
extern const char* exc_SwissCheese;
extern const char* exc_Limits;
extern const char* exc_Overflow;
extern const char* exc_HolesNotSet;
extern const char* exc_IDTooHigh; 
extern const char* exc_TooMany;
extern const char* exc_Infinite;

#define PLUS_OR_MINUS_ONE_TO_SYMMETRISE +1

/** chain-agnostic base for base **/
class BaseData {
public:
	// number strings of type
	vector<int> its_structure;
	// number of down spins (i.e. finite roots plus infinite raps.
	int number_magnons;						
	// number of holes (see definition below)
	/// FIXME: number_holes should be removed, it's only used by Base, but is somehow entangled in numberFreedoms... fix later.
	int number_holes;		
	int number_spinons;
public:	
	// create empty
	BaseData (void);
	// set everything
	BaseData (const int number_down, const vector<int>& structure, const int number_holes, const int number_spinons);
	// copy
	BaseData (const BaseData& original);
	// destroy
	~BaseData(void) {};
	// assign
	BaseData& operator= (const BaseData& rhs);
	
	/// NOTE: we mustn't check number_holes! (it sorta doesn't exist)
	inline bool operator== (const BaseData& another) const 
	{ return ( (number_magnons==another.number_magnons) && (number_spinons==another.number_spinons) && ( its_structure == another.its_structure ) ); };
	// two bases are unequal if they're not equal	
	inline bool operator!= (const BaseData& another) const { return !(*this==another); };
	
	
	// the basis for sorting bases.
	int numberFreedoms (void) const;
protected:	
	void crop(void);
};

/** base: number strings of every type, and the whole shebang **/

class Base: public BaseData {
public:
	// pointer to chain on which we live (make sure it doesn't get deleted!)
	const Chain* p_chain;
protected:							
	// allowed state ids, per sector, for this base 0 <= id < maxId . last element [number_sectors] is state maxid.
	vector<long long int> lim_id;			 
public:
	// number of slots for each sector. see below for the definition of a sector.
	vector<int> number_slots_sector;		
	// number of filled slots per sector // TODO: this member is almost an exact copy of the parent object. do something about it.
	vector<int> number_particles_sector; 	
	
public:	
	Base (const Chain* const on_chain, const BaseData& base_data);
	// create a base on the basis of a given structure
	Base (const Chain* const on_chain, const vector<int>& structure, const int number_holes=0);
	// create a base on the basis of a given structure and number of down spins
	Base (const Chain* const on_chain, const int number_down, const vector<int>& structure, const int number_holes=0);
	// create a base with only particles in the 1+interval and no holes (ground state base)
	Base (const Chain* const on_chain, const int the_number_roots);		
	// copy constructor
	Base (const Base& original); 
	// destructor
	virtual ~Base(void) {};
	// assign
	Base& operator= (const Base& rhs);
	// clone
	virtual Base* clone(void) const = 0;	
	
	// two bases are equal if they have equal chains, numbers down, numbers holes, and elements
	inline bool operator== (const Base& another) const
	{ return ( (*p_chain == *another.p_chain) && (number_magnons==another.number_magnons) && (number_holes==another.number_holes) && ( its_structure == another.structure() ) ); };
	// two bases are unequal if they're not equal	
	inline bool operator!= (const Base& another) const { return !(*this==another); };
	
	// member access
	Base& setStructure (const vector<int>& structure, const int new_number_spinons);
	Base& setNumberHoles (const int new_number_holes);
	Base& setNumberSpinons (const int number_spinons);
		
	vector<int> structure(void) const { return its_structure; };
	
	// get the id number of a set of shifts
	long long int id (const vector<Young>& shifts) const;
	// get the shifts belonging to a certain id number
	vector<Young> shifts(const long long int id) const;
	
	// by default, numberDown and numberRoots are the same.
	// however, if there are infinite rapidities, there are included in numberDown but not in numberRoots
	inline int numberDown (void) const { return number_magnons; };	
	virtual inline int numberRoots (void) const { return number_magnons; }; 
	inline int numberInfiniteRapidities (void) const { return numberDown() - numberRoots(); };

	// number of types with at least one instance (note that the constructor should make sure there are no trailing zeroes)
	inline int numberTypes (void) const { return its_structure.size(); };
	
	// what we are all about... the number of instances of each string type
	inline int numberStringsOfType (int string_type) const { return (its_structure.size()>string_type)?its_structure[string_type]:0; };
	
	// the lowest string type consists of two sectors, the others of one (see below)
	inline int numberSectors (void) const { return ((1==its_structure.size())&&(!number_holes))?1:its_structure.size()+1; };
	
	// number of slots for each sector. see below for the definition of a sector.
	inline vector<int> numberSlotsPerSector (void) const { return number_slots_sector; };
	
	// number of filled slots per sector
	inline vector<int> numberParticlesPerSector (void) const { return number_particles_sector; };
	
	// the number of slots per string type (depends on parity of string type as well as limQuantumNumbers)
	vector<int> numberSlotsPerStringType (void) const; 
	
	// last element: the lowest integer that is not an allowed id number
	// other elements describe the size of lower sectors (young tableaux) analogously
	inline vector<long long int> limId (void) const { return lim_id; };	
	
	
	// the lowest integer that is not an allowed quantum number:   ceil(2*I_max)   
	virtual vector<int> limQuantumNumbers (void) const = 0; 
	
	// a hole is an unfilled slot in the ground sector.
	// should the ground sector be larger than the number of particles of string type 1+, 
	// as can happen in the case of parity mismatch between that number and the number of slots for string type 1+,
	// (depending on the constant PLUS_OR_MINUS_ONE_TO_SYMMETRISE) there is always a hole in the ground sector. 
	inline int numberHoles (void) const { return number_holes; };
	
	// a spinon is a filled slot in sector[1]
	inline int numberSpinons (void) const { return number_particles_sector[1]; };
	
	// the minimum number of holes is the difference between the number of ground slots and 1+strings.
	inline int minNumberHoles (void) const { return max (0, number_slots_sector[0] - its_structure[0]); };
	inline int minNumberSpinons (void) const { return max (0, its_structure[0] - number_slots_sector[0]); };
		
	// the maximum number of holes: either the number of ground slots (duh) or the number of places the particles can go to.
	inline int maxNumberHoles (void) const { return min (number_slots_sector[0], number_slots_sector[1] - minNumberSpinons() ); }; 
	inline int maxNumberSpinons (void) const { return min (number_slots_sector[0] - minNumberHoles() , number_slots_sector[1] ); };
	
	
	
protected:
	void checkConsistency (void) ;
	void calculateNumberSlots(void);
	void calculateNumberParticles(void);
 	void calculateLimId (void);
};


		// ** What is a sector ? **
		//
		// every string type produces an REAL max_quantum_number. 
		// For string_type 1+, the (half-)integers I such that abs(I)< max_quantum_number constitute the (1+)-interval.
		// 	(1+)-interval[0] .. (1+)-interval[number_1+_interval - 1] 
		// For 1-, this produces the (1-)-interval, &c.
		
		// sector[0] is the part of the (1+)-interval that would be filled 
		// if the system were in the subspace ground state for the given number of down spins.
		// if the leftmost (most negative I) position in the (1+)-interval is (1+)-interval[0],
		// the leftmost position in sector[0], 		sector[0][0] == (1+)-interval[ground_sector_offset]
		// the rightmost position in sector[0] 		sector[0][number_slots[0]-1] == (1+)-interval[ground_sector_offset + number_slots[0] -1] 
		
		// sector[1] is what is left of the (1+)-interval. 
		// its leftmost position 			sector[1][0]== (1+)-interval[0]
		// going on up to and including 	sector[1][ground_sector_offset-1] == (1+)-interval[ground_sector_offset-1] 
		// skipping sector[0], 
		// continuing from and including 	sector[1][ground_sector_offset] == (1+)-interval[ground_sector_offset + number_slots[0]] 
		// up to and including				sector[1][number_StringsOfType(1) - number_slots[0] - 1] 
		
		// sector[2] equals the (1-)-interval, sector[3] the next string_type interval, &c.
		
		// a sector consists of slots, indicating the number of allowed quantum number values in that sector
		// it is filled with particles (a particle being a 1+string, 1-string, 2string, etc.)
		

		// note that the number of slots in the ground sector may be less than the number of down spins,
		// if the quantum number limits do not allow for so many slots
		// in that case, however, the number of holes is still the number of down spins minus 
		// the number of particles in the ground sector
		// some holes, therefore, are not represented by empty slots in the ground sector. 
		// we call these Dark Holes (tm).


#endif
