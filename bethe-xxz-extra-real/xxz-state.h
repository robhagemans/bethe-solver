#ifndef XXZ_STATE_H
#define XXZ_STATE_H

#include "bethe.h"
#include "strip.h"
#include "state.h"
#include "det.h"

extern const char* exc_NotXXZ;
 
	/** Gapless State **/
	class XXZ_State : public ::State {
	public:
		// set to subspace vacuum
		XXZ_State (const Base* const ground_base); 				
		// set to given id for given base (negative id: do not set id)
		XXZ_State (const Base* const on_base, const long long int state_id):  State(on_base, state_id), tfh_rapidity(on_base) { setTables(); setFreeRapidities(); };
		// set quantum numbers directly	
		XXZ_State (const Base* const on_base, const Strip<int>& quantum):  State(on_base, quantum), tfh_rapidity(on_base) { setTables(); setFreeRapidities(); };	
		// set shifts
		XXZ_State (const Base* const on_base, const vector<Young>& shift): 	State(on_base, shift), tfh_rapidity(on_base) { setTables(); setFreeRapidities(); }; 
		// copy constructor
		XXZ_State (const XXZ_State& original): 	State(original), tfh_rapidity(original.tfh_rapidity) { setTables(); };	
		// copy from generic state	
		XXZ_State (const State& original);
		// destructor
		virtual ~XXZ_State () { delete[] tf_half_n_anis; delete[] sf_half_n_anis; };
		// assignment
		//operator=(const State& original);
		//XXZ_State& operator= (const XXZ_State& rhs);
		
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
		virtual REAL plusFormFactor (const State& state_mu) const { throw Exception("XXZ_State::plusFormFactor", exc_NotImplemented); }
		
		// deviation from string hypothesis
		virtual REAL stringDeviation (void) const;	
		// Gaudin's matrix. sets norm.
		virtual Square<REAL> matrixGaudin (void) const;
		
		// log of determinant of H-2P matrix
		virtual REAL lndetHm2P (const State& state_mu) const;
		// log of determinant of H- matrix
		virtual REAL lndetHminus (const State& state_mu) const;
		
	protected:		
		// tanh of rapidity
		Strip<REAL> tfh_rapidity;	
		// sin of n\zeta/2
		REAL* sf_half_n_anis;
		// tan of n\zeta/2
		REAL* tf_half_n_anis;
		
		// set the optimisation tables
		void setTables(void);

		// calculate the energy and set the energy field
		REAL calculateEnergy (void) const;
		
		// set rapidities and related data fields 
		void setRapidity (const int j, const int alpha, const long double rapidity);
		void setRapidityAtIndex (const int index, const long double rapidity);
	private:	
		// internal methods for readability
		long double scatteringDerivative (const int j, const int alpha, const int k, const int beta) const;
		long double thetaDerivative (const long double rapidity, const int length, const int parity) const;
		long double lhsDerivative(const int j, const int alpha) const;
		
		// internal methods for readability
		long double rhs (const int j, const int alpha) const;
		long double thetaTfh(const long double th_rapidity, const int length, const int parity) const; 
		long double tfhInvTheta (const long double phase, const int string_type) const ;
	};

#endif
