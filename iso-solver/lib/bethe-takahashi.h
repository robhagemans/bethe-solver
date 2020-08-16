#ifndef BETHE_TAKAHASHI_H
#define BETHE_TAKAHASHI_H


#include <math.h>
#include <vector>
#include <complex>
#include <iostream>


using namespace std;

#include "basics.h"
#include "exception.h"
#include "strip.h"
#include "square.h"
#include "configuration.h"

/* exceptions */

extern const string exc_Runaway;
extern const string exc_NonFinite;


/** solver ABC **/
class BetheSolver {
	// solve for roots
	virtual bool solve() = 0;
	// all roots
	virtual vector< complex<double> > roots (void) const = 0;
	// one root
	virtual complex<double> root (const int j, const int alpha, const int a) const = 0;
};



/** Isotropic Bethe Solver **/

class IsoBetheTakahashiSolver : public BetheSolver {
public:
	/* general solution policies */

	// maximum iterations
	static int max_iterations;
	// required precision for convergence
	static double precision;


	/* newton's method policies */

	// maximum Newton's method steps
	static int max_newton;
	// precision to extrapolate to before newton's method is invoked
	static double newton_threshold;
	// factor with which to increase newton threshold if newton doesn't work
	static double newton_factor;
	// factor of worsening we allow from newton's method before stopping it
	static double newton_bandwidth;

	/* solve() strategy */

	// number of consecutive newton steps allowed
	static int newton_consecutive;
	// number of consecutive iterations/extrapolations allowed
	static int extrapolate_consecutive;

	/* policies for iterate() */

	// maximum number of steps when running
	static int run_max;
	// between 0 and 1; higher number means slower run.
	static float run_sloth;

	/* thresholds, epsilons, sentinels, etc. */

    // to check if we have real part zero in quantum number
	static double threshold_quantum_number;
	// indicates no convergence: some big number.
    static double no_convergence;


public:

	/* set data */

	// constructor
	IsoBetheTakahashiSolver(const int chainLength, const IsoConfiguration& config, const Strip<int>& the_quantum_numbers);
	// set quantum numbers directly
	void setQuantum2I (const Strip<int>& the_quantum_numbers);

 	// destructor
	virtual ~IsoBetheTakahashiSolver() {};


	/* the solver */

	// reset rapidities, convergence
	virtual void reset();
	// solve for rapidities and deviance
	virtual bool solve();
	// iterate rapidities
	virtual void iterateBT();


	/* readout */

	// roots of the bethe equation
	virtual vector< complex<double> > roots() const;
	// one root
	virtual complex<double> root (const int j, const int alpha, const int a) const;


	/* accuracy measures */

	// deviation from string hypothesis
 	virtual double stringDeviation() const;

	// error in solution of (product) Bethe equation: accuracy check.
	double betheError(void) const;


	/* classification */

	// check admissiblity (non-singularity)
	bool admissible() const;
	bool symmetric() const;

	bool singular() const;
	bool oddSymmetric() const;


	/* quantum numbers */

	// bethe quantum numbers (2I) from takahashi quantum numbers, only for two-strings but no lolution required
	vector< int > cleanQuantum2JFrom2I(void) const;
	// bethe J's from roots, unrounded
	vector< complex<double> > dirtyQuantumJ(void) const;

	static vector< complex<double> > dirtyQuantumJforRoots (const int for_chain_length, const vector< complex<double> >& for_roots);

	// bethe quantum numbers (2J) from roots
	vector< int > cleanQuantum2J (void) const;


protected:

    // iterate a few times, the extrapolate
    bool solveExtrapolate (const int max_iter, const double precision);

	// apply Newton's method once
	bool newton(double& newt_convergence);

	// Bethe equation -- used by newton()
	long double betheZero (const int j, const int alpha) const;
	// gaudin matrix - used by newton()
	Square<double> matrixGaudin (void) const;
	// scattering derivative - used by matrixGaudin
	long double scatteringDerivative (const int j, const int alpha, const int k, const int beta) const;


	// rhs of bethe equation -- used by iterate()
	long double rhs (const int j, const int alpha) const;
	// scattering term -- used by iterate()
	long double theta (const long double th_rapidity, const int length) const;




protected:
	// N, chain length
	int chain_length_;
	// string configuration
	IsoConfiguration config_;

	// Bethe-Takahashi quantum_number == 2*I
	Strip<int> quantum_number;
	// finite string rapidities
	Strip<double> rapidity;

	// square diff of last iteration
	double convergence;
	// number of Bethe-Takahashi iterations undergone
	int iterations;
	// number of Newton iterations undergone
	int newton_iterations;



};


#endif
