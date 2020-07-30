#ifndef ISO_DEVIATED_H
#define ISO_DEVIATED_H


#include <math.h>
#include <vector>
#include <complex>
#include <iostream>
using namespace std;

#include "strip.h"
#include "iso-state.h"



/** exceptions **/

extern const char* exc_NoCanDo;
extern const char* exc_WrongPairWidth;
extern const char* exc_Deprecated;



/** XXXDeviatedState **/

class XXXDeviatedState : public XXX_State {
public:
	// stores imaginary deviation
	FullStrip<REAL> deviance; 
	// stores real deviation
	FullStrip<REAL> aberration;
	
public:
	// maximum number of iteratons in deviance solution
	static int max_iter_deviance;   
	// precision to be reached in deviance iteration
	static REAL deviance_precision; 
	// how bad 'convergence' may become before we call it quits
	static REAL solve_bailout;
	// damping for deviance iteration
	static REAL damping_deviance;
	// damping for solveSymmetric()
	static REAL damping_symmetric;
	
	// deviance as set in constructor. can't be zero, that would give non-finite results.
	const static REAL start_deviance;
	// epsilon in the check for deviance=-0.5. NOTE that this need not be very small, as we don't really get 'near-collapses'	
	const static REAL epsilon_collapse;
	// below which the norm would give a non-finite result if we would use the complex roots (i.e. we should use the string equations for the norm)
	const static REAL epsilon_deviation_norm;
	// how close roots may come to i/2 before we consider it to be exactly that for the calculation of energies
	const static REAL epsilon_energy;

	
protected:
	// hold this value in iteration (rapidity)
	Strip<short int> hold;  
	FullStrip<short int> hold_deviation;  
	
public:	
	// we can currently only construct an XXXDeviatedState from a State 
	XXXDeviatedState (const State& original);
	virtual ~XXXDeviatedState () {};
	
	// iterate rapidities, taking deviance into account
	virtual void iterate (void);
	// solve for rapidities and deviance
	virtual bool solve(void);		
	
	// norm of state 
	virtual REAL norm(void) const;
	// roots of the bethe equation, with deviance for twostrings
	virtual vector< complex<REAL> > roots (void) const;
	// one root
	virtual complex<REAL> root (const int j, const int alpha, const int a) const;
	
	/// Gaudin's matrix -- used by norm() and newton() 
//	Square<complex<REAL> > matrixGaudin (void);

	// we cannot calculate form factors with a deviated state as right state (lambda); use the converse form factor
	virtual REAL longitudinalFormFactor (const State& state_mu) const;
	virtual REAL transverseFormFactor (const State& state_mu) const; 
	virtual REAL plusFormFactor (const State& state_mu) const;
	
	// deviation from string hypothesis
 	virtual REAL stringDeviation (void) const;		
 	
 	// magnitude of deviance
	virtual REAL devianceMagnitude(void) const;
	
	// is this an odd symmetric state?
	bool oddSymmetric(void) const;
protected:		
	// calculate energy and set energy field
	virtual REAL calculateEnergy (void) const;
	
	// steps in iteration: get an iteration of the rapidity
	virtual long double newRapidity (const int j, const int alpha) const;
	// steps in iteration: get a value for deviance and iteration
	virtual bool getNewDeviances (FullStrip<REAL>& new_deviance, FullStrip<REAL>& new_aberration, const int j, const int alpha) const;
	
	// solve some classes of symmetric states
	bool solveSymmetric(void);
	bool solveSingular(void);
	
	// these are needed by norm()
	// d small theta/d lambda 
	virtual complex<long double> thetaDerivative (const complex<long double> rap, const int length) const;
	// scattering derivative term (d big theta / d lambda) 
	virtual complex<long double> scatteringDerivative (const int j, const int alpha, const int a, const int k, const int beta, const int b) const;
	virtual complex<long double> scatteringDerivativeLeftDev (const int j, const int alpha, const int a, const int length_k, const int k, const int beta) const;
	virtual complex<long double> scatteringDerivativeRightDev (const int undev_j, const int j, const int alpha, const int k, const int beta, const int b) const;
	virtual long double scatteringDerivativeNoDev (const int undev_j, const int j, const int alpha, const int undev_k, const int k, const int beta) const;
	
	
// RLH2008			
	double findRoot1(const int bethe_2xj, const double eps, const double lam_guess);
	double exRealIndicator1(const int bethe_2xj, const double eps, const double lam);
	double exRealIndicator2(const int bethe_2xj, const double eps, const double lam);
	bool findExtraReal(const int j, const int alpha, const int a);
	
	int getSum2xBetheQuantum (const int j, const int alpha) const;
	void rebase(const int j, const int alpha, const int a);
	
};


#endif
