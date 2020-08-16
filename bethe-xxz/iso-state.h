#ifndef ISO_STATE_H
#define ISO_STATE_H

#include "bethe.h"
#include "state.h"
extern const char* exc_NotXXX;


	/** Isotropic State **/
	class XXX_State : public ::State {
	public:
		// set to subspace vacuum
		XXX_State (const Base* const ground_base);
		// set to given id for given base
		XXX_State (const Base* const on_base, const long long int state_id): State(on_base, state_id) { setFreeRapidities(); };
		// set quantum numbers directly
		XXX_State (const Base* const on_base, const Strip<int>& quantum):  State(on_base, quantum) { setFreeRapidities(); };
		// set shifts
		XXX_State (const Base* const on_base, const vector<Young>& shift): 	State(on_base, shift) { setFreeRapidities(); };
		// copy constructor
		XXX_State (const XXX_State& original): State(original) {};
		// copy from generic state
		XXX_State (const State& original): State(original)	{ if (original.p_chain->delta() != 1.0)	throw Exception("XXX_State::XXX_State", exc_NotXXX); };
		// destructor
		virtual ~XXX_State () {};
		// assignment
		//virtual XXX_State& operator= (onst XXX_State& rhs);
		//operator=(const State& original);

		// get finite complex roots
		virtual vector< complex<REAL> > roots(void) const;
		// get one finite complex root
		inline virtual complex<REAL> root (const int j, const int alpha, const int a) const
		{	return	rapidity(j, alpha) + 0.5 * ( p_chain->stringLength(j) - 1 - 2*a ) *I;	};

		// reset rapidities, convergence
		virtual void setFreeRapidities(void);
		// Bethe equation
		virtual long double betheZero (const int j, const int alpha) const;
		// iterate Bethe equations once
		virtual void iterate (void);

		// form factor <state_mu | Sz | *this> ; sets norm field if uncalculated
		virtual REAL longitudinalFormFactor(const State& state_mu) const;
		// form factor <state_mu | S- | *this> ; sets norm field if uncalculated
		virtual REAL transverseFormFactor(const State& state_mu) const;
		// form factor <state_mu | S+ | *this> ; sets norm field if uncalculated
		virtual REAL plusFormFactor(const State& state_mu) const;

		// deviation from string hypothesis
		virtual REAL stringDeviation(void) const;

		// reduced Gaudin's matrix. sets norm.
		virtual Square<REAL> matrixGaudin(void) const;
		// log of determinant of H-2P matrix
		virtual REAL lndetHm2P(const State& state_mu) const;
		// log of determinant of H- matrix
		virtual REAL lndetHminus(const State& state_mu) const;

		///the naming of these three methods is horrible
		// bethe I's from rapidities
		vector< complex<REAL> > calculateBetheI(void) const;
		// bethe quantum numbers (2I) from rapidities
		vector< int > calculateBethe2I (void) const;
		// bethe quantum numbers (2I) from takahashi quantum numbers, only for two-strings!
		vector< int > calculateBetheQuantumNumbers(void) const;
		// error in solution of (product) Bethe equation: accuracy check.
		REAL XXX_State::betheError(void) const;

	protected:
		// set rapidities and related data fields
		inline virtual void setRapidity (const int j, const int alpha, const long double new_rapidity)
		{ 	rapidity(j, alpha) = new_rapidity; };
		inline virtual void setRapidityAtIndex (const int index, const long double new_rapidity)
		{	rapidity.element(index) = new_rapidity;	};

		// calculate the energy and set the energy field
		virtual REAL calculateEnergy (void) const;

		// internal methods for readability
		inline long double thetaDerivative (const long double rap, const int length) const
		{	return 4.0 / ( 1.0*length  + 4.0*rap*rap/(1.0*length) ); 	};	// also used by iso-deviated
		// internal methods for readability
		long double scatteringDerivative (const int j, const int alpha, const int k, const int beta) const;// also used by iso-deviated
	private:
		long double rhs (const int j, const int alpha) const;
		long double theta (const long double th_rapidity, const int length) const;
	};

#endif
