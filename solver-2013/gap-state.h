#ifndef GAP_STATE_H
#define GAP_STATE_H

#include "bethe.h"
#include "strip.h"
#include "state.h"
#include "det.h"

extern const char* exc_RapidityMinHalfPi;
extern const char* exc_RapidityCollision;

	/** Gapped State **/
	class Gap_State : public ::State {
	public:
		// set to subspace vacuum
		Gap_State (const Base* const ground_base);
		// set to given id for given base
		Gap_State (const Base* const on_base, const long long int state_id);
		// set quantum numbers directly
		Gap_State (const Base* const on_base, const Strip<int>& quantum);
		// set shifts
		Gap_State (const Base* const on_base, const vector<Young>& shift);
		// copy constructor
		Gap_State (const Gap_State& original): State(original), tfh_rapidity(original.tfh_rapidity), ising_rapidity_over_pi(original.ising_rapidity_over_pi)
		{	setTables(); 	};
		// destructor
		virtual ~Gap_State () {  delete[] tf_half_n_anis; delete[] sf_half_n_anis;  };
		//operator=(const State& original);
		//Gap_State& operator= (const Gap_State& rhs);

		// set id. returns distance (square diff) with previous state
		virtual int setId(const long long int state_id);
		// check admissiblity
		virtual bool admissible(void) const;

		// get finite complex roots
		virtual complex<REAL> root (const int j, const int alpha, const int a) const;
		// get one finite complex root
		virtual vector< complex<REAL> > roots (void) const;

		// reset rapidities, convergence
		virtual void setFreeRapidities(void);
		// Bethe equation
		virtual long double betheZero (const int j, const int alpha) const;
		// iterate Bethe equations once
		virtual void iterate (void);

		// form factor <state_mu | Sz | *this> ; sets norm field if uncalculated
		virtual REAL longitudinalFormFactor (const State& state_mu) const;
		// form factor <state_mu | S- | *this> ; sets norm field if uncalculated
		virtual REAL transverseFormFactor (const State& state_mu) const;
		// form factor <state_mu | S+ | *this> ; sets norm field if uncalculated
		virtual REAL plusFormFactor (const State& state_mu) const { throw Exception("Gap_State::plusFormFactor", exc_NotImplemented); }

		// deviation from string hypothesis
		virtual REAL stringDeviation (void) const;

		// reduced Gaudin's matrix. sets norm.
		virtual Square<REAL> matrixGaudin (void) const;
		// log of determinant of H-2P matrix
		virtual REAL lndetHm2P (const State& state_mu) const;
		// log of determinant of H- matrix
		virtual REAL lndetHminus (const State& state_mu) const;


	protected:
		// tan of rapidity
		Strip<REAL> tfh_rapidity;
		// rapidity in the Ising limit
		Strip<REAL> ising_rapidity_over_pi;

		// sinh of n\zeta/2
		REAL* sf_half_n_anis;
		// tanh of n\zeta/2
		REAL* tf_half_n_anis;


		// set the optimisation tables
		void setTables(void);
		void setIsingRapidities (void);

		// set rapidities and related data fields
		void setRapidity (const int j, const int alpha, const long double rapidity);
		void setRapidityAtIndex (const int index, const long double rapidity);

		// calculate the energy and set the energy field (called by energy())
		REAL calculateEnergy (void) const;
	private:
		// internal methods for readability
		long double scatteringDerivative (const int j, const int alpha, const int k, const int beta) const;
		long double thetaDerivative (const long double lambda, const int length) const;

		// internal methods for readability
		long double rhs (const int j, const int alpha) const;
		long double thetaTfh(const long double the_tfh_rapidity, const REAL the_rapidity, const int length) const;
		long double invTheta (const long double phase, const int string_type) const ;
	};

#endif
