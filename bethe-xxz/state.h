#ifndef STATE_H
#define STATE_H

#include "bethe.h"
#include "strip.h"
#include "square.h"
#include "det.h"
#include "chain.h"
#include "base.h"
#include "young.h"


/** exceptions **/

extern const char* exc_NotGroundBase;
extern const char* exc_NegativeHoles;
//extern const char* exc_Forbidden; //State::SetID
extern const char* exc_Strings;
extern const char* exc_NumberDown;
extern const char* exc_Runaway;
extern const char* exc_NonFinite;
extern const char* exc_Uncontrol;
extern const char* exc_Equal;	// matrixHm2P(), matrixHminus()
extern const char* exc_NotImplemented;


/**  sentinel values **/

// initial value for convergence; some big number.
extern const REAL NO_CONVERGENCE;
// sentinel value to flag uncalculated REAL fields
extern const REAL NOT_CALCULATED;

// don't set state id
#define NO_ID -1
// no mode
#define NO_MODE -1
// integer field not set
#define NOT_SET -666



class State  {
public:
	const Chain* p_chain; 					// the chain on which this is a state
	const Base* p_base;						// the base of this state

//protected:
	Strip<int> quantum_number;				// quantum_number == 2*I
	Strip<REAL> rapidity;					// finite string rapidities

	REAL convergence;						// square diff of last iteration
	int iterations;							// number of Bethe iterations undergone
	int newton_iterations;					// number of Newton iterations undergone

protected:
	// to be set by functions requesting them, to be set to NOT_CALCULATED if rapidities change:
	mutable REAL its_lnnorm;				// stores norm
	mutable REAL its_energy;				// stores energy
	///mutable REAL its_deviation;	///TODO
	///mutable short int admissible;

	// to be set by functions requesting them, to be set to NOT_CALCULATED if quantum numbers change:
	mutable long long int its_id;			// stores id
	mutable vector<Young> its_shifts;		// stores shifts
	mutable int its_mode;					// stores mode number (== momentum index)
	mutable short int its_symmetry;			// stores symmetry.

public:
	// general solution policies
	static int max_iterations;				// maximum iterations
	static REAL precision;					// required precision for convergence
	// newton's method policies
	static int max_newton;					// maximum Newton's method steps
	static REAL newton_threshold;			// precision to extrapolate to before newton's method is invoked
	static REAL newton_factor;				// factor with which to increase newton threshold if newton doesn't work
	static REAL newton_bandwidth;			// factor of worsening we allow from newton's method before stopping it
	// solve() strategy
	static int newton_consecutive;			// number of consecutive newton steps allowed
	static int extrapolate_consecutive;		// number of consecutive iterations/extrapolations allowed
	// policies for iterate()
	static int run_max;						// maximum number of steps when running
	static float run_sloth; 				// between 0 and 1; higher number means slower run.


public:
	// set to given id for given base (default: 'ground state')
	State (const Base* const on_base, const long long int state_id=0);
	// copy constructor
	State (const State& original);
	// set shifts
	State (const Base* const on_base, const vector< Young >& shift);
	// set quantum numbers directly
	State (const Base* const on_base, const Strip<int>& quantum);
	// destructor
	virtual ~State () {};

	// assignment
	//virtual State& operator= (const State& rhs);

	// set id. returns distance (square diff) with previous state
	inline virtual int setId(const long long int state_id) 	{ 	return setShifts (p_base->shifts(its_id=state_id)); 	};
	// set shifts. returns distance (square diff) with previous state
	int setShifts(const vector<Young>& shift);
	// set single shifts
	void setShift(const int sector, const Young& shift);
	// set quantum numbers directly
	void ::State::setQuantumNumbers (const Strip<int> the_quantum_numbers);
	// check admissiblity
	virtual bool admissible(void) const;

	// get/calculate the ID of this state
	inline long long int id(void) const { return (its_id==NO_ID)? its_id = p_base->id(shifts()) : its_id; };
 	// get the particle/hole shifts for this state
	inline vector<Young> shifts(void) const { return its_shifts.size()? its_shifts : calculateShifts();  };
 	// returns whether the state is symmetric
 	inline bool symmetric(void) const { return (its_symmetry==NOT_SET)? its_symmetry = calculateSymmetry() : its_symmetry; };

	// mode number, between 0 (inclusive) and chain.length() (exclusive)
	inline int mode(void) const { return (NO_MODE==its_mode)? calculateMode() : its_mode; };
	// momentum of the state, between 0 (inclusive) and 2 pi (exclusive)
	inline REAL momentum(void) const { return PI*mode() / (0.5*p_chain->length()); };
	// energy (calculated first time, then stored)
	inline REAL energy(void) const { return (its_energy==NOT_CALCULATED)? calculateEnergy(): its_energy ; };



	// reset rapidities, convergence
	virtual void setFreeRapidities(void) =0;
	// iterate Bethe equations once
	virtual void iterate(void) = 0;
	// iterate Newton's method once
	void newton(REAL& newt_convergence);
	// plain iteration to convergence
	virtual bool solve(const int max_iter, const REAL precision);
	// iteration and interpolation to convergece
	virtual bool solveExtrapolate(const int max_iter, const REAL precision);
	// iteration, extrapolation, newton's method to convergence
	virtual bool solve(void);

	// number of holes in the Fermi interval
	inline int numberHoles(void) const { return p_base->numberHoles(); };
	// number of (finite) string rapidities
	inline int numberRapidities(void) const { return quantum_number.numberElements(); };

	// get finite complex roots
	virtual vector< complex<REAL> > roots(void) const = 0;
	// get one finite complex root
	virtual complex<REAL> root(const int j, const int alpha, const int a) const =0;

	// Gaudin's matrix. must set its_lnnorm -- used by norm() and newton()
	virtual Square<REAL> matrixGaudin(void) const =0;
	// norm of the state (calculated the first time, then stored) // matrixGaudin sets lnnorm. discard its result.
	inline virtual REAL norm(void) const { return (  ((its_lnnorm==NOT_CALCULATED) ? matrixGaudin() :0),  its_lnnorm  );	};


	// first-order deviation from string hypothesis
	virtual REAL stringDeviation(void) const = 0;
	// magnitude of deviation in exactly-deviated result (zero for B--T solutions)
	virtual inline REAL devianceMagnitude(void) const { return 0.0; };

	// equals zero if Bethe's equation satisfied -- used by newton()
	virtual long double betheZero(const int j, const int alpha) const = 0;

	// form factor <state_mu | Sz | *this> ; sets norm field if uncalculated
	virtual REAL longitudinalFormFactor(const State& state_mu) const = 0;
	// form factor <state_mu | S- | *this> ; sets norm field if uncalculated
	virtual REAL transverseFormFactor(const State& state_mu) const = 0;
	// form factor <state_mu | S+ | *this> ; sets norm field if uncalculated
	virtual REAL plusFormFactor(const State& state_mu) const = 0;

	// log of determinant of H-2P matrix
	virtual REAL lndetHm2P(const State& state_mu) const = 0;
	// log of determinant of H- matrix
	virtual REAL lndetHminus(const State& state_mu) const = 0;

protected:
	// calculate number of holes
	int calculateNumberHoles(void) const;
	// calculate shifts and set mutable shifts field
	vector<Young> calculateShifts(void) const;
	// calculate mode and set mutable mode field
	int calculateMode(void) const;
	// calculate energy and set energy field
	virtual REAL calculateEnergy (void) const = 0;
	// calculate symmetry and set symmetry field
	bool calculateSymmetry(void) const;

	// set a rapidity and dependent fields (such as tanh storage)
	virtual void setRapidity (const int j, const int alpha, const long double rapidity) = 0;
	virtual void setRapidityAtIndex (const int index, const long double rapidity) = 0;

};

#endif
