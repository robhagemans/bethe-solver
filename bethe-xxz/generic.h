#ifndef GENERIC_H
#define GENERIC_H

#include "bethe.h"
#include "matrix.h"

#include "xxz-chain.h"
#include "iso-chain.h"
#include "gap-chain.h"

#include "xxz-base.h"
#include "iso-base.h"
#include "gap-base.h"

#include "xxz-state.h"
#include "iso-state.h"
#include "gap-state.h"

#include "iso-deviated.h"


// generally: this value is not set
#define NOT_SET -666


/** exceptions **/
extern const char* exc_Ferro;				// generic contructors
extern const char* exc_InfiniteRaps;
extern const char* exc_NegativeAnisotropy;	// name()
extern const char* exc_BadName;
extern const char* exc_ThresholdOne;
extern const char* exc_NullPointer;
extern const char* exc_NotConverged; 		//solve()
extern const char* exc_Deviated; 			// solve()


/** policies **/

// default cutoff for number of string types allowed
#define CUTOFF_TYPES 4


/** file names **/
// a file name consists of moniker-value pairs started by a Quantity.name(), separated by sep_Name and terminated by sep_Final plus an extension

// separators
extern const string sep_Name;
extern const string sep_Final;
extern const string sep_Vector;

// monikers
extern const string mon_File;	// file format moniker
extern const string mon_Base;
extern const string mon_Delta;
extern const string mon_Length;
extern const string mon_LeftDown;
extern const string mon_RightDown;
extern const string mon_Spinon;

// number of characters needed to display positive integer number
inline int fieldWidth (const int number) { return 1+int( floor(log(1.0*number)/log(10.0)) ); };

// convert monikers and values into names
string nameFromValue (const string moniker, const REAL value);
string nameFromValue (const string moniker, const int value, const int width=0);
string nameFromValue (const string moniker, const vector<int>& value);

// find the value corresponding to the given moniker in a string (which may contain more moniker-value pairs)
int intValueFromName(const string name, const string moniker);
REAL realValueFromName(const string name, const string moniker);
string valueFromName(const string name, const string moniker);

// the name of a base
string name(const BaseData& base_data);
string name(const Base* const p_base);
// the name of a chain. if (left) number_down != 0, a magnetisation name is added
string name(const Chain* const p_chain, const int number_down=0);

/** generic construction **/

// copy chain
inline Chain* copy (const Chain* const p_chain) { return p_chain->clone(); };
// copy base
inline Base* copy (const Base* const p_base) { return p_base->clone(); };
// copy state
State* copy (State* const p_state) ;

// convert a name into a new chain
Chain* newChain (const string name, const int cutoff_types=CUTOFF_TYPES);
// convert a name into a new base (and a new chain if p_chain equals zero)
Base* newBase (const string name, Chain*& p_chain);
BaseData readBaseData (const string name);

// create a new chain form given data
Chain* newChain (const double delta, const int chain_length, const int cutoff_types = CUTOFF_TYPES);
// create a new base
Base* newBase (const Chain* const p_chain, const BaseData& base_data);
Base* newBase (const Chain* const p_chain, const vector<int>& base_vec, const int number_spinons, const int number_infinite_rapidities = 0);
Base* newBase (const Chain* const p_chain, const int number_down, const vector<int>& base_vec, const int number_spinons, const int number_infinite_rapidities = 0);
// create a ground base with a given number of down spins
Base* newGroundBase (const Chain* const chain, const int number_down);
// create a ground state for a given ground base
State* newGroundState (const Base* const base);
// create a state with a given id. if id==-1, the id is not set.
State* newState (const Base* const base, const long long int id) ;
// create a state with given shifts
State* newState (const Base* const base, const vector<Young>& shifts) ;


// create all bases ordered by particle content and create a list of pointers to them
vector<BaseData> allNewBases (
	const Chain* const p_chain, const int number_down,
	const int max_string_length,			// cutoff string length
	const int uptoinc_number_particles, 	// cutoff number particles
	const int uptoinc_number_spinons, 		// cutoff number spinons
	const int max_infinite = 2
);


// generic solve
State* solve(State* const p_state, const REAL deviation_threshold, const bool deviate = false);

#endif
