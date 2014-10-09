#ifndef ISO_SOLVER_H
#define ISO_SOLVER_H

#include <math.h>
#include <vector>
#include <complex>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

#include "basics.h"
#include "exception.h"
#include "strip.h"
#include "bethe-takahashi.h"

/** exceptions **/

extern const string exc_Runaway;
extern const string exc_NonFinite;



class IsoSolver : public IsoBetheTakahashiSolver {
	
public: 

	// construct
	IsoSolver(const int chainLength, const IsoConfiguration& config, const Strip<int>& the_quantum_numbers);
	
	// destruct
	virtual ~IsoSolver() {};
	
	/* solver */
	
	// reset rapidities, convergence
	virtual void reset(); 
	// solve for rapidities and deviance
	virtual bool solve();		
	// iterate rapidities, taking deviance into account
	virtual void iterate();
	
	/* roots readout */
	
	// roots of the bethe equation
	virtual vector< complex<double> > roots (void) const;
	// one root
	virtual complex<double> root (const int j, const int alpha, const int a) const;
	
	/* accuracy measures */
	
	// magnitude of deviance
	double devianceMagnitude() const;
	
	// deviation from string hypothesis
 	virtual double stringDeviation() const;	
 	
	/* classification */
	
	
	// sum of bethe's quantum numbers (from takahashi qn)
    int sumOfBethe2xJ(const int j, const int alpha);

	
public: 	
 	/* solver policies */
 	
 	// maximum string deviation we accept before trying to solve for deviance
	static double deviation_threshold;
	
	// maximum number of iteratons in deviance solution
	static int max_iter_deviance;   
	// precision to be reached in deviance iteration
	static double deviance_precision; 
	// how bad 'convergence' may become before we call it quits
	static double solve_bailout;
	// damping for deviance iteration
	static double damping_deviance;
	// damping for solveSymmetric()
	static double damping_symmetric;
	
	/* sentinels */

	// signals absence of convergence
	const static double no_convergence;
	// signals failure to calculate a double-valued function
	const static double not_calculated;
	
protected:	
	/* internal settings */
	
	// deviance as set in constructor. can't be zero, that would give non-finite results.
	const static double start_deviance;
	// epsilon in the check for deviance=-0.5. NOTE that this need not be very small, as we don't really get 'near-collapses'	
	const static double epsilon_collapse;
	// below which the norm would give a non-finite result if we would use the complex roots (i.e. we should use the string equations for the norm)
	const static double epsilon_deviation_norm;
	// how close roots may come to i/2 before we consider it to be exactly that for the calculation of energies
	const static double epsilon_energy;

	
	
protected:
    void resetDeviances();

 	// plain iteration to convergence 
	bool solveIterate(const int max_iter, const double precision); 
	// iteration and interpolation to convergece
	bool solveExtrapolate(const int max_iter, const double precision); 

	// solve some classes of symmetric states
	bool solveSymmetric();
	bool solveSingular();

	// steps in iteration: get an iteration of the rapidity
	long double newRapidity (const int j, const int alpha) const;
	// steps in iteration: get a value for deviance and iteration
	bool getNewDeviances (FullStrip<double>& new_deviance, FullStrip<double>& new_aberration, const int j, const int alpha) const;

protected:
// RLH2008			
	double findRoot1(const int bethe_2xj, const double eps, const double lam_guess);
	double exRealIndicator1(const int bethe_2xj, const double eps, const double lam);
	double exRealIndicator2(const int bethe_2xj, const double eps, const double lam);
	bool findExtraReal(const int j, const int alpha, const int a);


protected:
	// stores imaginary deviation
	FullStrip<double> deviance; 
	// stores real deviation
	FullStrip<double> aberration;

	// hold this value in iteration (rapidity)
	Strip<short int> hold;  
	FullStrip<short int> hold_deviation;  

};





#endif
