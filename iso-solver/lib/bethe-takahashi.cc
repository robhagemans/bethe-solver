#include "bethe-takahashi.h"


/* default solution policy (works best for short chains) */

/* general */
// maximum iterations
int IsoBetheTakahashiSolver::max_iterations	= 1000;		
// required precision for convergence
double IsoBetheTakahashiSolver::precision = 1e-28;	
/* newton's method */
// maximum Newton's method steps
int IsoBetheTakahashiSolver::max_newton = 100;		
// precision to extrapolate to before newton's method is invoked
double IsoBetheTakahashiSolver::newton_threshold = 1e-2;	
// factor with which to increase newton threshold if newton doesn't work	
double IsoBetheTakahashiSolver::newton_factor	= 1e-2;		
// factor of worsening we allow from newton's method before stopping it
double IsoBetheTakahashiSolver::newton_bandwidth = 10.0;		
/* solve() strategy */
// number of consecutive newton steps allowed
int IsoBetheTakahashiSolver::newton_consecutive	= 7;
// number of consecutive iterations/extrapolations allowed	
int IsoBetheTakahashiSolver::extrapolate_consecutive  = 20;	
/* policies for iterate()  */
// maximum number of steps when running
int IsoBetheTakahashiSolver::run_max = 50;		
// between 0 and 1; higher number means slower run.
float IsoBetheTakahashiSolver::run_sloth = 0.5; 		


/* thresholds, epsilons, sentinels, etc. */

// to check if we have real part zero in quantum number
double IsoBetheTakahashiSolver::threshold_quantum_number=1e-8;
// indicates no convergence: some big number.
double IsoBetheTakahashiSolver::no_convergence = sq(BOOMBANOOMBA)*4.0; 


/* exception descriptions thrown in this module */
const string exc_Runaway = "runaway rapidity";
const string exc_NonFinite = "non-finite result";



// defines sorting order for quantum numbers (for use with std::sort() )
// true if a is less than b according to the ordering:
// 0, -2, 2, -4, 4, -6, 6, etc.
// -1, 1, -3, 3, -5, 5, etc.
inline bool quantumNumberLessThan (const int a, const int b) 
{ 	return (abs(a)<abs(b)) || ( (abs(a)==abs(b)) && (a<b) ); 	}





/** interpolate **/
// source: Numerical Recipes
 // Polynomial interpolation/extrapolation, NR page 113.
void polint (vector<double>& xa, vector<double>& ya, const double x, double& y, double& dy)
{
  int i, m, ns=0;
  double den, dif, dift, ho, hp, w;

  int n = xa.size();
  vector<double> c(n), d(n);
  dif = fabs(x-xa[0]);
  for (i = 0; i < n; i++) {
    if ((dift = fabs(x-xa[i])) < dif) {
      ns = i;
      dif = dift;
    }
    c[i] = ya[i];
    d[i] = ya[i];
  }
  y = ya[ns--];
  for (m = 1; m < n; m++) {
    for (i = 0; i < n-m; i++) {
      ho = xa[i] - x;
      hp = xa[i+m] - x;
      w = c[i+1] - d[i];
      if ((den = ho-hp) == 0.0) throw Exception ("polint", "zero denominator");
      den = w/den;
      d[i] = hp * den;
      c[i] = ho * den;
    }
    y += (dy = (2 * (ns + 1) < (n-m) ? c[ns + 1] : d[ns--]));
  }
}


    IsoBetheTakahashiSolver::IsoBetheTakahashiSolver (const int chain_length, const IsoConfiguration& on_base, const Strip<int>& quantum)
        :   chain_length_(chain_length), config_(on_base), quantum_number(quantum),
		    rapidity(on_base), convergence(no_convergence), iterations(0), newton_iterations(0)
    { 	
	    reset(); 
    }




    /* set quantum numbers */
    void IsoBetheTakahashiSolver::setQuantum2I (const Strip<int>& the_quantum_numbers)
    {
	    const char* here = "::State::setQuantumNumbers";

	    // copy the quantum numbers
	    quantum_number = the_quantum_numbers;
	
	    // indicate uncalculated fields
	    convergence = IsoBetheTakahashiSolver::no_convergence;	
	    iterations = 0;
	    newton_iterations = 0;
    }

	
	
	
	/* set free rapidies */
	void IsoBetheTakahashiSolver::reset()
	{
		// simply the bethe equation without scattering phases
		// j iterates over string types; alpha iterates over all instantiations of a string type; i.e. over all strings.
		for (int j=0; j < config_.numberTypes(); ++j) 
		for (int alpha=0; alpha < config_.numberStringsOfType(j); ++alpha) 
			rapidity(j, alpha) = 0.5*config_.stringLength(j) * tan(0.5*PI*quantum_number(j, alpha)/(1.0*chain_length_));
cerr<<quantum_number<<endl;
cerr<<rapidity<<endl;		
		iterations=0;
		newton_iterations=0;
		convergence = IsoBetheTakahashiSolver::no_convergence; 
	}


	
	
	
	
	/* roots of the bethe equation */
	vector< complex<double> > IsoBetheTakahashiSolver::roots (void) const
	{
		vector< complex<double> > root_vec (config_.numberRoots());
		for (int index=0, j=0; j < config_.numberTypes(); ++j)  
		for (int alpha=0; alpha < config_.numberStringsOfType(j); ++alpha) 
		for (int a=0; a < config_.stringLength(j); ++a, ++index) {
			root_vec[index] =  root(j, alpha, a); 
		}
		return root_vec;
	}


    /* roots of the bethe equation */
    complex<double> IsoBetheTakahashiSolver::root (const int j, const int alpha, const int a) const
    {	
	    return rapidity(j, alpha) 	+  ( 0.5 *(config_.stringLength(j)-1-2*a) )*I;	
    }


	/* one iteration of Bethe-Takahashi equation */
	void IsoBetheTakahashiSolver::iterateBT (void) 
	{
		string here = "IsoBetheTakahashiSolver::iterate"; 

		Strip<double> new_rapidity (config_);
		double square_diff = 0.0;	
		double one_over_n = 1.0/(1.0*chain_length_);
		long double new_rap = 0.0;
		
		for (int j=0; j < config_.numberTypes(); ++j) 
		for (int alpha=0; alpha < config_.numberStringsOfType(j); ++alpha) {
			double current_rapidity = rapidity(j, alpha);	// this should be finite & <1
			for (int run = 0; run < run_max; ++run) {
				long double last_rap = new_rap;
				new_rap = 0.5*config_.stringLength(j) * tan( 0.5*one_over_n*IsoBetheTakahashiSolver::rhs (j, alpha) );
				if (finite(new_rap)) break;
				// if it doesn't work, slowly raise the offending rapidity and see if we can get something reasonable.
				// NOTE: since rhs only uses tfh_rapidity, we don't calculate the atan. 
				// thus rapidity[][] is temporarily out of sync !!
				rapidity(j, alpha) *= 1.1;
				if (!finite(rapidity(j, alpha))) break; // give up.
			}
			rapidity(j, alpha) = current_rapidity; // reset the old value, don't mess with current rapidities.
			new_rapidity(j, alpha) = new_rap;
			if (!finite(new_rapidity(j, alpha))) throw Exception (here, exc_Runaway);
			square_diff += sq(rapidity(j, alpha) - new_rapidity(j, alpha));
		}
		// if we have NaN square diff, throw exception. *this shouldn't happen.*
		if  (!finite(square_diff)) throw Exception (here, exc_NonFinite); 
				
		// set the new values for the rapidities
		for (int j=0; j < config_.numberTypes(); ++j) 
		for (int alpha=0; alpha < config_.numberStringsOfType(j); ++alpha) { 
			rapidity(j, alpha) = new_rapidity(j, alpha);
		}
		convergence = square_diff; 
		++iterations;
	}



	/* deviation from string hypothesis: sum of squares of norm of deltas */
	double IsoBetheTakahashiSolver::stringDeviation (void) const
	{
		double sum_deviations=0;
		for (int index=0, j=0; j < config_.numberTypes(); ++j)  {
			int length_j = config_.stringLength(j);
			if (1==length_j) continue; // no contribution from 1-strings: no deltas.
			
			for (int alpha=0; alpha < config_.numberStringsOfType(j); ++alpha) {
				complex<double> deviation [length_j];
				for (int a=0; a < length_j-1; ++a, ++index) 
				{
					complex<double> root_j = root(j, alpha, a);
					deviation[a] = pow(  (root_j + 0.5*I) / (root_j - 0.5*I), - chain_length_);

					// product over other strings
					for (int k=0; k < config_.numberTypes(); ++k) 
					for (int beta=0; beta < config_.numberStringsOfType(k); ++beta) {
						if ((k==j) && (alpha==beta)) continue;
						for (int b=0; b< config_.stringLength(k); ++b) {
							complex<double> root_k = root (k, beta, b); 
							deviation[a] *= (root_j - root_k + I) / (root_j - root_k - I);
						}
					}
					
					// product over the same string
					int s1 = 2*(length_j-a);
					int s2 = 2*(length_j-a-1);
					int s3 = 2*a;
					int s4 = 2*(a-1);
					deviation[a] *= 1.0*(s1?s1:1)*(s2?s2:1)/( (s3?s3:1)*(s4?s4:1) );

					if (a>0) deviation[a] *= deviation[a-1];
					sum_deviations += std::norm(deviation[a]);
				}
			}
		}
		return sqrt(sum_deviations);
	}
	
	
	
	/* calculate the Bethe quantum numbers (2J), from J = I1+I2+N/2, I1-I2 = [wide] */
	vector< int > IsoBetheTakahashiSolver::cleanQuantum2JFrom2I (void) const
	{
		const int twostring_type=1;
		
		// currently, we ignore infinite rapidities
		vector< int > bethe_quantum_number (config_.numberRoots());
		for (int index_j_alpha = 0, j=0; j< config_.numberTypes(); ++j)
		for (int alpha=0; alpha < config_.numberStringsOfType(j); ++alpha) 
		for (int a=0; a< (j+1) ; ++a, ++index_j_alpha) {
			if (0==j) bethe_quantum_number[index_j_alpha] = quantum_number(j, alpha);
			else if (twostring_type==j) {
				// criterion for narrow/wide pair, from I$ == I1 + I2 + N/2 ; I1 - I2 = [wide] 
				bool narrow = ( (-chain_length_ - config_.numberStringsOfType(twostring_type)+1 + quantum_number(twostring_type, alpha))/2 + config_.numberRoots() )%2;
				bethe_quantum_number[index_j_alpha] = quantum_number(j, alpha)/2 + chain_length_/2 + (narrow?0:(a?1:-1));	
			}
			else throw Exception ("IsoBetheTakahashiSolver::cleanQuantum2JFrom2I", "not implemented for strings longer than 2");
		}
		return bethe_quantum_number;
	}
	
		

	/* calculate the quantum numbers (J) as they are in the complex (non-Takahashi) Bethe equations, from the solution */
	vector< complex<double> > IsoBetheTakahashiSolver::dirtyQuantumJ (void) const
	{
/*
		// currently, we ignore infinite rapidities
		vector< complex<double> > bethe_roots = roots();
		vector< complex<double> > bethe_i (config_.numberRoots());
		for (int alpha=0; alpha < config_.numberRoots(); ++alpha) {
			if (abs(real(bethe_roots[alpha])) < threshold_quantum_number && abs(abs(imag(bethe_roots[alpha]))-0.5)< threshold_quantum_number ) {
				// numberRoots is number of finite roots. if even, there is no zero; if odd, there is a zero.
				int root_sign = isgn(imag(bethe_roots[alpha]));
				if (config_.numberRoots()%2) 
					// odd. there's a zero q.n.
					bethe_i[alpha] = 0.25*(chain_length_ - root_sign*(chain_length_%4));
				else 
					// even. no zero.
					bethe_i[alpha] = 0.25*(chain_length_ - root_sign*(2 - chain_length_%4));
				continue;
			}
			complex<double> lhs = atan(2.0*bethe_roots[alpha]) * double(chain_length_);
			for (int beta=0; beta<config_.numberRoots(); ++beta){ 
				if (alpha!=beta) lhs -= atan(bethe_roots[alpha] - bethe_roots[beta]);
			}
			bethe_i[alpha] = lhs/PI;	
		}
		return bethe_i;
*/
        return dirtyQuantumJforRoots(chain_length_, roots());
	}
	
	
	/* calculate the quantum numbers (J) as they are in the complex (non-Takahashi) Bethe equations, from the solution */
	vector< complex<double> > IsoBetheTakahashiSolver::dirtyQuantumJforRoots (const int for_chain_length, const vector< complex<double> >& for_roots)
	{
	 	vector< complex<double> > bethe_i (for_roots.size());
		for (int alpha=0; alpha < for_roots.size(); ++alpha) {
			if (abs(real(for_roots[alpha])) < threshold_quantum_number && abs(abs(imag(for_roots[alpha]))-0.5)< threshold_quantum_number ) {
				// numberRoots is number of finite roots. if even, there is no zero; if odd, there is a zero.
				int root_sign = isgn(imag(for_roots[alpha]));
				if (for_roots.size()%2) 
					// odd. there's a zero q.n.
					bethe_i[alpha] = 0.25*(for_chain_length - root_sign*(for_chain_length%4));
				else 
					// even. no zero.
					bethe_i[alpha] = 0.25*(for_chain_length - root_sign*(2 - for_chain_length%4));
				continue;
			}
			complex<double> lhs = atan(2.0*for_roots[alpha]) * double(for_chain_length);
			for (int beta=0; beta<for_roots.size(); ++beta){ 
				if (alpha!=beta) lhs -= atan(for_roots[alpha] - for_roots[beta]);
			}
			bethe_i[alpha] = lhs/PI;	
		}
		return bethe_i;
	}
	
	/* calculate the integer quantum numbers (2J) as they are in the complex (non-Takahashi) Bethe equations, from the solution. round to int */
	vector< int > IsoBetheTakahashiSolver::cleanQuantum2J (void) const
	{
		vector< complex<double> > bethe_i = dirtyQuantumJ();
		vector<int> bethe_quantum_number (bethe_i.size());
		for (int i=0; i<bethe_i.size();++i) 
			bethe_quantum_number[i] = int( round(2.0*real(bethe_i[i])) );
		return bethe_quantum_number;
	}	
	
	
	
	/* check the complex Bethe equations */
	double IsoBetheTakahashiSolver::betheError(void) const
	{
		vector<complex<double> > bethe_j = dirtyQuantumJ();
		double error = 0.0;
		for (int i=0; i<bethe_j.size();++i) 
			error += ::norm(bethe_j[i] -  0.5*round(2.0*real(bethe_j[i])));
		return error;

	}


	
	
	
	/* sum of scattering terms */
	inline long double IsoBetheTakahashiSolver::theta (const long double rap, const int length) const 
	{	return 2.0 * atan( 2.0 * rap / (1.0*length) ) ; 	}


	/* bethe equation: right hand side */
	long double IsoBetheTakahashiSolver::rhs (const int j, const int alpha) const
	{
		long double scattering_term = 0.0;
		int length_j = config_.stringLength(j);
		
		for (int k=0; k < config_.numberTypes(); ++k) 
		for (int beta=0; beta < config_.numberStringsOfType(k); ++beta) {
			int length_k = config_.stringLength(k);
			// epsilon of long double equals that of intermediate results
			long double lj_m_lk = rapidity(j, alpha) - rapidity(k, beta);
			scattering_term += IsoBetheTakahashiSolver::theta (lj_m_lk, length_j + length_k);
			if (length_j != length_k) 
				scattering_term += IsoBetheTakahashiSolver::theta (lj_m_lk, abs(length_j-length_k));
			for (int i=1; i< min(length_j, length_k); ++i) 
				scattering_term += 2.0 * IsoBetheTakahashiSolver::theta (lj_m_lk, abs(length_j-length_k) + 2*i);
		}
		return PI*quantum_number(j, alpha) + scattering_term;
	}

	/* Bethe equation */
	inline long double IsoBetheTakahashiSolver::betheZero (const int j, const int alpha) const
	{	
		return chain_length_ * IsoBetheTakahashiSolver::theta (rapidity(j, alpha), config_.stringLength(j)) - IsoBetheTakahashiSolver::rhs(j, alpha);	
	};
		
		
/* these are originally from class State; they are Delta-agnostic */
		
		
    /** check admissibility of state **/
    bool IsoBetheTakahashiSolver::admissible() const
    {
	    return !singular();
    }
	
	
    bool IsoBetheTakahashiSolver::singular() const
    {
	    // check quantum numbers. a few states need to be excluded 
	    // because the string hypothesis leads to non-finite results
	
	    // * the quantum numbers are symmetrically configured and there exist even-string-length bound states.
	    //   - there will be an even string with double rapidity zero;
	    //   - this implies (on the LHS of the original bethe equations) a lambda_j that is +- I zeta/2,
	    //     yielding a zero numerator or denominator on the left, where no infinities should occur.
	    //     (only (1+something)^N -type divergencies, being cancelled out by lambda_j - lambda_k = i zeta + delta on the right.)
	    //  Solutions in this class are actually proper Bethe states (given the correct limiting procedure, see notes). 
	    //  However, their contribution to correlations is necessarily zero (see notes).
	
	    // * we have symmetrically distributed quantum numbers and there is more than one type of odd string of the same parity 
	    //   e.g. a three-string with positive parity as well as the regular one-plusses.
	    //   - we can have more than one zero rapidity (doublely, only at zero). this is forbidden by the exclusion principle.
	    // 	NOTE: actually, this latter class is solvable by deviation. Also, its solutions do contribute to correlations.
	    // 	So it must only be excluded if we don't deviate. Therefore, these should be counted as 'deviated', not 'inadmissible'
			
	    // check for a zero in an even-length state. assume ordered quantum numbers.
	    for (int j=0; j < config_.numberTypes(); ++j) {
		    if (config_.numberStringsOfType(j) && !quantum_number(j, 0) && !(config_.stringLength(j)%2)) 
			    return symmetric();
	    }	
	    return false;
    }
	

	/** find 'odd symmetric' states which cause problemas but have nonzero contributions **/
	bool IsoBetheTakahashiSolver::oddSymmetric() const
	{
		//  no symmetry no problem
		if (!symmetric()) return false;
		// assumes ordered quantum numbers.
		// count the number of odd-length strings with quantum number zero
		int number_odd_zero = 0;
		for (int j=0; j<config_.numberTypes(); ++j) 
			if (config_.numberStringsOfType(j)>0 && !quantum_number(j,0) && (config_.stringLength(j)%2)) ++number_odd_zero;			
		
		// more than one odd zero: problem.
		return (number_odd_zero>1);
	}
		

    bool IsoBetheTakahashiSolver::symmetric(void) const
    {
	    // let's see if the distribution of quantum numbers is symmetric.
	    /*
	    // this should work for unordered quantum numbers
	    for (int j=0; j < base.numberTypes(); ++j) { 
		    int sum_quantum_numbers = 0, sum_cube_quantum_numbers = 0;
		    for (int alpha=0; alpha < base.numberStringsOfType(j); ++alpha) {
			    sum_quantum_numbers += quantum_number(j, alpha);
			    sum_cube_quantum_numbers += quantum_number(j, alpha)*quantum_number(j, alpha)*quantum_number(j, alpha);
		    }
		    if (sum_quantum_numbers || sum_cube_quantum_numbers) return false;
	    }
	    return true;
	    */
	
	    // faster: assume ordering 0 -1 1 -2 2 etc which we enforce anyway
	    for (int j=0; j < config_.numberTypes(); ++j) { 
		    int number_j = config_.numberStringsOfType(j);
		    // central string not zero: not symmetric
		    if (number_j%2 && quantum_number(j,0)) return false;
		    for (int alpha=number_j%2; alpha < number_j; alpha+=2) {
			    // if not symmetric, bail out immediately
			    if (quantum_number(j, alpha+1) != -quantum_number(j, alpha)) return false;
		    }
	    }
	    return true;
    }		

    /** one iteration of Newton's method **/
    bool IsoBetheTakahashiSolver::newton (double& newt_convergence)
    {
	    // this is probably not generic...
	    const double tfh_threshold = 10.0;	
	
	    // Newton's method: Numerical Recipes in C,  p381
	    int sign_permutation = 1;
	    vector<int> indx (rapidity.numberElements()); // used by backsub, given by decompose 
	    vector<double> delta (rapidity.numberElements());
	    newt_convergence = 0.0;
	
	    // Gaudin's matrix is Jacobi's matrix for Bethe's equations 
	    Square<double> jacobian = matrixGaudin(); 

	    // delta is the right side of   -F in mx_J - vc_deltax == - vec_F
	    for (int j=0; j < config_.numberTypes(); ++j) 
	    for (int alpha=0; alpha < config_.numberStringsOfType(j); ++alpha) 
		    delta[rapidity.index(j,alpha)] = betheZero(j, alpha);

	    if (!decomposeLU (jacobian, indx, sign_permutation)) return false;
	    backsubLU (jacobian, indx, delta);
	
	    for (int i=0; i < rapidity.numberElements(); ++i) {
		    if (abs(delta[i]) > tfh_threshold) continue; // don't generate tanhs that are 1
		    if (isNan(delta[i])) continue;
		    rapidity.element(i) = rapidity.element(i) + delta[i]; 
		    newt_convergence += sq(delta[i]);
	    }
	
	    ++newton_iterations;
	    // iterate once to check convergence
	    iterateBT();	
	    return true;		
    }
	
	/*** ***/
	
	// internal methods for readability
	inline long double thetaDerivative (const long double rap, const int length) 
	{	return 4.0 / ( 1.0*length  + 4.0*rap*rap/(1.0*length) ); 	};	
	
	/* scattering derivative term (d big theta / d lambda) */	
	long double IsoBetheTakahashiSolver::scatteringDerivative (const int j, const int alpha, const int k, const int beta) const
	{
		int length_j = config_.stringLength(j);
		int length_k = config_.stringLength(k);
		//int prod_parity = p_chain->string_parity[j]*p_chain->string_parity[k];
		double lj_m_lk = rapidity(j, alpha) - rapidity(k, beta);
		long double big_theta = thetaDerivative (lj_m_lk, length_j+length_k);
		if (length_j != length_k) 
			big_theta += thetaDerivative (lj_m_lk, abs(length_j - length_k));
		for (int i=1; i< min(length_j, length_k); ++i) 
			big_theta += 2.0 * thetaDerivative (lj_m_lk, abs(length_j - length_k) + 2*i);
		
		return big_theta;
	}
	
	
	/* Gaudin's matrix */
	Square<double> IsoBetheTakahashiSolver::matrixGaudin (void) const
	{
		const char* here = "IsoBetheTakahasiSolver::matrixGaudin";
		Square<double> gaudin (rapidity.numberElements()); 
		
		for (int j=0; j< config_.numberTypes(); ++j)
		for (int alpha=0; alpha < config_.numberStringsOfType(j); ++alpha) {
		
			int index_j_alpha = rapidity.index(j, alpha);	
			for (int k=0; k< config_.numberTypes(); ++k)
			for (int beta=0; beta < config_.numberStringsOfType(k); ++beta) {
				
				int index_k_beta = rapidity.index(k, beta);		
				if (index_j_alpha == index_k_beta) {
					gaudin[index_j_alpha][index_k_beta] = 
						chain_length_ * thetaDerivative(rapidity(j, alpha), config_.stringLength(j)); 
	
					for (int l=0; l< config_.numberTypes(); ++l)
					for (int gamma=0; gamma < config_.numberStringsOfType(l); ++gamma) {
						if (rapidity.index(l, gamma) != index_j_alpha)
							gaudin[index_j_alpha][index_k_beta] -= scatteringDerivative (j, alpha, l, gamma);
					}
				}
				else 	gaudin[index_j_alpha][index_k_beta] = scatteringDerivative (j, alpha, k, beta);
				if (!finite(gaudin[index_j_alpha][index_k_beta])) throw Exception (here, exc_NonFinite);
			}
		}
		
		// set the norm value now that we're at it
//		its_lnnorm = lndet(gaudin);
		return gaudin;
	}
	
	
/** simply iterate until convergence is reached **/
/*
bool IsoBetheTakahashiSolver::solveIterate (const int max_iter, const double precision) 
{
	while ( (convergence > precision) && (++iterations < max_iter) ) iterate(); 
	return (convergence <= precision);
}
*/


/** iterate a few times, then extrapolate  **/
bool IsoBetheTakahashiSolver::solveExtrapolate (const int max_iter, const double precision) 
{
	const double rapidity_threshold = 18.0; 
	const int number_extrapolate = 4;  
	double dummy_convergence = 0.0;
	double last_convergence = no_convergence;
	Strip<double>* rapidity_store[number_extrapolate];
	for (int i=0; i<number_extrapolate; ++i) {
		rapidity_store[i] = new Strip<double>(config_);
	}
	vector<double> xvalues (number_extrapolate);
	vector<double> yvalues (number_extrapolate);
	while ( (convergence > precision) && (iterations < max_iter) ) {
cerr<<iterations<<" "<<rapidity<<endl;		
		for (int i=0; i<number_extrapolate; ++i){
			iterateBT();
			if (convergence <= precision) break;
			*(rapidity_store[i]) = rapidity; 		//copy rapidity strip
		}
		if (convergence <= precision) break;
		last_convergence=convergence;
		// extrapolate along a trial-and-error-deduced magic curve
		vector<double> yvalues (number_extrapolate);
		for (int i=0; i< number_extrapolate; ++i) 
		    xvalues[i] = 2.0/(1.0*(i+4));	// the magic curve
		
		for (int j=0; j < rapidity.numberElements(); ++j) {
			for (int i=0; i<number_extrapolate; ++i) yvalues[i] = (*rapidity_store[i]).element(j);
			polint (xvalues, yvalues, 0.0, rapidity.element(j), dummy_convergence);

			// don't extrapolate to unreasonable values, fall back to last iteration.
			if (!finite(rapidity.element(j)) || (abs(rapidity.element(j)) > rapidity_threshold)) 
				rapidity.element(j) = (*rapidity_store[number_extrapolate-1]).element(j);
		}
		// check convergence
		iterateBT (); 
		// if extrapolation didn't help..
		if (convergence > last_convergence) 
			for (int j=0; j < rapidity.numberElements(); ++j) 
			    rapidity.element(j) = (*rapidity_store[number_extrapolate-1]).element(j) ;
	}
	for (int i=0; i<number_extrapolate; ++i) {
		 delete rapidity_store[i];
	}
	return (convergence <= precision);
}


/** iterate, extrapolate, newton's method: the whole shebang **/
bool IsoBetheTakahashiSolver::solve ()
{
	double newt_convergence = no_convergence;
	int last_iter_newton = 0;
	// first solve up to some precision
	solveExtrapolate (max_iterations, newton_threshold);
	while ( (convergence > precision) && (iterations < max_iterations)  ) {
		double last_convergence = convergence;
		// try newton's method within a given fluctuation band 
		// for at most a given number of times
		double last_newt_convergence = 2.0*newt_convergence;  
		last_iter_newton = newton_iterations;
		while ((newt_convergence < newton_bandwidth * last_newt_convergence)
				&& newton_iterations < max_newton
				&& newton_iterations < last_iter_newton + newton_consecutive
				) {
			last_newt_convergence = newt_convergence;
			// try a newton step
			if (!newton (newt_convergence)) 
			{
			    // didn't work.
				// must be decomposeLU: singular matrix.
				// newton can't be used, see if iterating once helps
//				its_lnnorm = NOT_CALCULATED;
				iterateBT ();
			}
			
			if (convergence<=precision) break; 
//			else its_lnnorm = NOT_CALCULATED;
		}
		if (convergence<=precision) break;
		
		// if we're here, newton failed for now.
//		its_lnnorm = NOT_CALCULATED;
		
		if (newton_iterations >= max_newton) {
			// we have reached the max. no other choice but to try interpolation
			// if we reach max_iter, we haven't converged.
			solveExtrapolate (max_iterations, precision);
			break;
		}
		else {
			// we're here, so we haven't reached the max.
			// however, newton is not an option right now.

			// go a given factor from the current convergence.
			// do only a limited number of iterations, then newton again.
			newton_threshold = max(newton_factor * convergence, precision);
			solveExtrapolate (iterations + extrapolate_consecutive, newton_threshold);
		}
	}
	return (convergence <= precision);
}



