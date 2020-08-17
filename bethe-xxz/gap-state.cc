#include "gap-state.h"

#define EQUAL_THRESHOLD 1e-10
#define COLLISION_THRESHOLD 1e-10


const char* exc_RapidityMinHalfPi = "rapidity at -pi/2";
const char* exc_RapidityCollision = "rapidity collision";


	/** subspace ground state constructor **/
	Gap_State::Gap_State (const Base* const ground_base)
		: 	::State(ground_base), tfh_rapidity(ground_base), ising_rapidity_over_pi(ground_base)
	{
		//check consistency
		for (int i=1; i<p_base->numberTypes(); ++i)
			if (p_base->numberStringsOfType(i)) throw Exception ("Gap_State::Gap_State", exc_NotGroundBase);
		setTables();
		setIsingRapidities();
		setFreeRapidities();
	}


	/** given ID constructor **/
	Gap_State::Gap_State (const Base* const the_base, const long long int state_id)
		: 	::State(the_base, state_id), tfh_rapidity(the_base), ising_rapidity_over_pi(the_base)
	{
		setTables();
		setIsingRapidities();
		setFreeRapidities();
	}

	/** set shifts **/
	Gap_State::Gap_State (const Base* const the_base, const vector<Young>& shift)
	: 	::State(the_base, shift), tfh_rapidity(the_base), ising_rapidity_over_pi(the_base)
	{
		setTables();
		setIsingRapidities();
		setFreeRapidities();
	}

	/** set quantum numbers directly **/
	Gap_State::Gap_State (const Base* const on_base, const Strip<int>& quantum)
		: 	::State(on_base, quantum), tfh_rapidity(on_base), ising_rapidity_over_pi(on_base)
	{
		setTables();
		setIsingRapidities();
		setFreeRapidities();
	}



	/** set tables **/
	void Gap_State::setTables(void)
	{
		// last string is longest string
		int table_size = 2* p_chain->stringLength(p_chain->numberTypes() -1);
		tf_half_n_anis = new REAL[ table_size ];
		sf_half_n_anis = new REAL[ table_size ];
		for (int n=0; n < table_size; ++n) {
			tf_half_n_anis[n] = tanh (0.5*n*p_chain->zeta());
			sf_half_n_anis[n] = sinh (0.5*n*p_chain->zeta());
		}
	}

	/** set rapidities and related data fields **/
	/* from given rapidity */
	void Gap_State::setRapidity (const int j, const int alpha, const long double new_rapidity)
	{
		rapidity(j,alpha) = new_rapidity;
		tfh_rapidity(j,alpha) = tan(new_rapidity);
	}

	void Gap_State::setRapidityAtIndex (const int index, const long double new_rapidity)
	{
		rapidity.element(index) = new_rapidity;
		tfh_rapidity.element(index) = tan(new_rapidity);
	}

	/** initial state for iteration **/
	void Gap_State::setFreeRapidities (void)
	{
		// simply the bethe equation without scattering phases
		// j iterates over string types; alpha iterates over all instantiations of a string type; i.e. over all strings.

		for (int j=0; j < p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha) {
			rapidity(j,alpha) = Gap_State::invTheta ( PI*quantum_number(j,alpha)/(1.0*p_chain->length()), j);
			tfh_rapidity(j,alpha) = tan(rapidity(j,alpha));
			if (!finite(tfh_rapidity(j,alpha))) {
				rapidity(j,alpha) = 0.0;
				tfh_rapidity(j,alpha) = 0.0;
			}
		}
		iterations=0;
		newton_iterations=0;
		convergence = NO_CONVERGENCE;
	}


	/** set quantum numbers from an id & reset the ising rapidities **/
	int Gap_State::setId(const long long int state_id) {
		State::setId(state_id);
		setIsingRapidities();
	}

	/** initialise vector of rapidities in the Ising limit. to be called from constructors. **/
	void Gap_State::setIsingRapidities (void)
	{
		REAL total_sum_rapidities = 0.0;	// K tilde
		vector<REAL> sum_quantum (p_base->numberTypes());
		for (int type_j=0; type_j < p_base->numberTypes(); ++type_j) {
			for (int alpha=0; alpha < p_base->numberStringsOfType(type_j); ++alpha) {
				sum_quantum[type_j] += 0.5*quantum_number(type_j, alpha);
			}
			total_sum_rapidities += sum_quantum[type_j];
		}
		total_sum_rapidities /= p_chain->length();

		// determine K tilde sub j, the sum of rapidities of each separate string type.
		REAL running_sum = 0.0;
		vector<REAL> sum_rapidities (p_base->numberTypes()); // K tilde sub j
		for (int type_j=0; type_j < p_base->numberTypes()-1; ++type_j) {
			int n_j = p_chain->stringLength(type_j);
			REAL denominator = p_chain->length();
			for (int type_k = type_j-1; type_k >=0 ; --type_k) {
				int n_k = p_chain->stringLength(type_k);
				sum_rapidities[type_j] += (n_j - n_k) * sum_rapidities[type_k];
				// sum_rapidities[type_k] has already been determined since type_k<type_j
			}
			for (int type_k=0; type_k < p_base->numberTypes(); ++type_k) {
				denominator -= 2.0 * min(n_j,p_chain->stringLength(type_k)) * p_base->numberStringsOfType(type_k);
			}
			sum_rapidities[type_j] *= 2.0*p_base->numberStringsOfType(type_j);
			sum_rapidities[type_j] += sum_quantum[type_j] - 2.0*n_j*p_base->numberStringsOfType(type_j)*total_sum_rapidities;
			sum_rapidities[type_j] /= denominator;
			running_sum += sum_rapidities[type_j];
		}
		// highest string type: Ktilde j is what's left.
		sum_rapidities[p_base->numberTypes()-1] = total_sum_rapidities - running_sum;

		for (int type_j=0; type_j < p_base->numberTypes(); ++type_j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(type_j); ++alpha) {
			// calculate the rapidity (divided by PI) in the Ising limit
			int denominator = p_chain->length() + p_base->numberStringsOfType(type_j);
			for (int type_k=0; type_k < p_base->numberTypes(); ++type_k)
				denominator -= 2*p_base->numberStringsOfType(type_k)*min(p_chain->stringLength(type_j), p_chain->stringLength(type_k));
			REAL numerator = 0.5*quantum_number(type_j, alpha) - 2.0*p_chain->stringLength(type_j) * total_sum_rapidities + sum_rapidities[type_j];
			for (int type_k = type_j-1; type_k >=0 ; --type_k)
				numerator += 2.0 * (p_chain->stringLength(type_j) - p_chain->stringLength(type_k)) * sum_rapidities[type_k];
			ising_rapidity_over_pi (type_j, alpha) = numerator/(1.0*denominator);
		}
	}

	/** admissibility check for gapped states is quite different from XXZ gapless, XXX **/
	bool Gap_State::admissible(void) const
	{
		const char* here = "Gap_State::admissible";

		// select a single period
		for (int j=0; j < p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha)
			// accept everything bethween -pi/2 and pi/2, both boundaries inclusive
			// throw away -pi/2 later, only if it indeed converged to something below or at pi/2.
			if ((ising_rapidity_over_pi(j,alpha) - 0.5 > EQUAL_THRESHOLD) || (ising_rapidity_over_pi(j,alpha) + 0.5 < -EQUAL_THRESHOLD)) {
				return false;
			}


		// check for 'classic' inadmissible states
		// * even string with real part of rapidity zero (mod PI)
		//   - this implies (on the LHS of the original bethe equations) a lambda_j that is +- I zeta/2,
		//     yielding a zero numerator or denominator on the left, where no infinities should occur.
		//     (only (1+something)^N -type divergencies, being cancelled out by lambda_j - lambda_k = i zeta + delta on the right.)
		// * more than one zero root (real and imaginary parts) zero (mod PI); happens only in odd strings
		//	 - forbidden by the exclusion principle.
		int number_odd_zero = 0;
		for (int type_j=0; type_j < p_base->numberTypes(); ++type_j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(type_j); ++alpha) {
			// if this fraction is close to an integer, the rapidity is close to a multiple of pi.
			if (abs(rounderr (ising_rapidity_over_pi(type_j, alpha))) < COLLISION_THRESHOLD) {
				if (p_chain->stringLength(type_j)%2) {
					// count odd zeroes (mod PI)
					++number_odd_zero;
					if (number_odd_zero > 1) return false;
				}
				else return false;
			}
		}

		// we survived all checks!
		return true;

	}


	/** roots of the bethe equation **/
	vector< complex<REAL> > Gap_State::roots (void) const
	{
		vector< complex<REAL> > root (p_base->numberRoots());
		for (int index=0, j=0; j < p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha)
		for (int a=0; a < p_chain->stringLength(j); ++a, ++index) {
			root[index] =  rapidity(j,alpha)
				+ 0.5 * I * (p_chain->zeta() * (p_chain->stringLength(j) - 1 - 2*a) );
		}
		return root;
	}

	/** roots of the bethe equation **/
	inline complex<REAL> Gap_State::root (const int j, const int alpha, const int a) const
	{
		return	rapidity(j,alpha)	+ 0.5 * I * (p_chain->zeta() * (p_chain->stringLength(j) - 1 - 2*a) );
	}


	/** energy, in units of J **/
	REAL Gap_State::calculateEnergy (void) const
	{
		its_energy = 0.0;
		for (int string_type=0; string_type< p_base->numberTypes(); ++string_type)
			for (int instance=0; instance<p_base->numberStringsOfType(string_type); ++instance)
				its_energy +=	( sf_half_n_anis[2] * sf_half_n_anis[2*p_chain->stringLength(string_type)] )
								/
								(	cos( 2.0* rapidity(string_type, instance) )
									- cosh( p_chain->stringLength(string_type) *p_chain->zeta() )
								);
		return its_energy;
	}

	/** calculate Gaudin's matrix for this state and set the norm field **/
	/* kinetic derivative term (d theta / d lambda) */
	long double Gap_State::thetaDerivative(const long double lambda, const int length) const
	{
		return ( sf_half_n_anis[2*length]/(sq(sin(lambda)) + sq(sf_half_n_anis[length])) );
	}
	/* scattering derivative term (d big theta / d lambda) */
	long double Gap_State::scatteringDerivative (const int j, const int alpha, const int k, const int beta) const
	{
		int length_j = p_chain->string_length[j];
		int length_k = p_chain->string_length[k];
		double lj_m_lk = rapidity(j,alpha) - rapidity(k,beta);
		long double big_theta = Gap_State::thetaDerivative (lj_m_lk, length_j+length_k);
		if (length_j != length_k)
			big_theta += Gap_State::thetaDerivative (lj_m_lk, abs(length_j - length_k));
		for (int i=1; i< min(length_j, length_k); ++i)
			big_theta += 2.0 * Gap_State::thetaDerivative (lj_m_lk, abs(length_j - length_k) + 2*i);

		return big_theta;
	}
	/* Gaudin's matrix */
	Square<REAL> Gap_State::matrixGaudin (void) const
	{
		const char* here = "State::matrixGaudin";
		Square<REAL> gaudin (rapidity.numberElements());

		for (int j=0; j< p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha) {

			int index_j_alpha = rapidity.index(j, alpha);
			for (int k=0; k< p_base->numberTypes(); ++k)
			for (int beta=0; beta < p_base->numberStringsOfType(k); ++beta) {

				int index_k_beta = rapidity.index(k, beta);
				if (index_j_alpha == index_k_beta) {
					gaudin[index_j_alpha][index_k_beta] = 	p_chain->length() * Gap_State::thetaDerivative(rapidity(j, alpha), p_chain->stringLength(j));

					for (int l=0; l< p_base->numberTypes(); ++l)
					for (int gamma=0; gamma < p_base->numberStringsOfType(l); ++gamma) {
						if (rapidity.index(l, gamma) != index_j_alpha)
							gaudin[index_j_alpha][index_k_beta] -= Gap_State::scatteringDerivative (j, alpha, l, gamma);
					}
				}
				else 	gaudin[index_j_alpha][index_k_beta] = Gap_State::scatteringDerivative (j, alpha, k, beta);
				if (!finite(gaudin[index_j_alpha][index_k_beta])) throw Exception (here, exc_NonFinite);
			}
		}
		// set the norm value now that we're at it
		its_lnnorm = lndet(gaudin);
		return gaudin;
	}


	/** small theta **/
	inline long double Gap_State::thetaTfh(const long double the_tfh_rapidity, const REAL the_rapidity, const int length) const
	{
        /// shouldn't that be ceil (the_rap/PI -0.5) ?? We want to include PI/2 but exclude -PI/2...
		return 2.0 * atan( the_tfh_rapidity / tf_half_n_anis[length] ) + 2.0*PI*floor( the_rapidity/PI + 0.5);
	}

	/** rapidity from phase **/
 	inline long double Gap_State::invTheta (const long double phase, const int string_type) const
	{
		return  atan( tan(0.5*phase) * tf_half_n_anis[p_chain->string_length[string_type]] ) + PI*floor(0.5*phase/PI + 0.5);
	}

	/* this doesn't work because the information for the floor function (i.e. the choice of branch) is lost.
		inline long double Chain::tanRapidityFromTranslationPhase (long double phase, int string_type) const
	{	return  tan(0.5*phase) * tanh_half_n_anis[string_length[string_type]] ;	}
	*/

	/** bethe equation: right hand side **/
	long double Gap_State::rhs (const int j, const int alpha) const
	{
		long double scattering_term = 0.0;
		int length_j = p_chain->stringLength(j);
		int parity_j = p_chain->string_parity[j];
		REAL current_tan_rapidity = tfh_rapidity(j,alpha);

		for (int k=0; k < p_base->numberTypes(); ++k)
		for (int beta=0; beta < p_base->numberStringsOfType(k); ++beta) {
			int length_k = p_chain->stringLength(k);
			int prod_parity = parity_j*p_chain->string_parity[k];
			// epsilon of long double equals that of intermediate results
			long double t_lj_m_lk =
					(current_tan_rapidity - tfh_rapidity(k,beta))
					/ (1.0 + current_tan_rapidity*tfh_rapidity(k,beta));	/// + for tan !
			long double lj_m_lk = rapidity(j,alpha) - rapidity(k,beta);
			scattering_term += Gap_State::thetaTfh (t_lj_m_lk, lj_m_lk, length_j + length_k);
			if (length_j != length_k)
				scattering_term += Gap_State::thetaTfh (t_lj_m_lk, lj_m_lk,  abs(length_j-length_k));
			for (int i=1; i< min(length_j, length_k); ++i)
				scattering_term += 2.0 * Gap_State::thetaTfh (t_lj_m_lk, lj_m_lk, abs(length_j-length_k) + 2*i);
		}
		return PI*quantum_number(j,alpha) + scattering_term;
	}

	/** bethe equation **/
	inline long double Gap_State::betheZero (const int j, const int alpha) const
	{
		return p_chain->length() * Gap_State::thetaTfh (tfh_rapidity(j,alpha), rapidity(j,alpha), p_chain->stringLength(j)) - Gap_State::rhs (j, alpha);
	}

	/** one iteration of Bethe's equation **/
	void Gap_State::iterate (void)
	{
		string here = "Gap_State::iterate";
		Strip<REAL> new_rapidity (p_base);
		REAL square_diff = 0.0;
		REAL one_over_n = 1.0/(1.0*p_chain->length());
		for (int j=0; j < p_base->numberTypes(); ++j) {
			for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha) {
				new_rapidity(j,alpha) = Gap_State::invTheta (one_over_n * Gap_State::rhs (j, alpha), j);
				square_diff += sq(rapidity(j,alpha) - new_rapidity(j,alpha));
			}
		}
		// if we have NaN rapidities, throw exception. *this shouldn't happen.*
		if  (isNan(square_diff)) throw Exception (here, exc_NonFinite);
		// set the new values for the rapidities
		for (int j=0; j < p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha) {
			rapidity(j,alpha) = new_rapidity(j,alpha);
			tfh_rapidity(j,alpha) = tan(new_rapidity(j,alpha));
		}
		convergence = square_diff;
		++iterations;
	}


	/** deviation from string hypothesis: sum of squares of norms of deltas **/
	REAL Gap_State::stringDeviation (void) const
	{
		REAL sum_deviations=0;
		for (int index=0, j=0; j < p_base->numberTypes(); ++j)  {
			int length_j = p_chain->stringLength(j);
			if (1==length_j) continue; // no contribution from 1-strings: no deltas.

			for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha) {
				complex<REAL> deviation [length_j];
				for (int a=0; a < length_j-1; ++a, ++index)
				{
					complex<REAL> root_j = root(j, alpha, a);
					deviation[a] = pow( sin(root_j + 0.5*I*p_chain->zeta())/sin(root_j - 0.5*I*p_chain->zeta()), - p_chain->length());

					// product over other strings
					for (int k=0; k < p_base->numberTypes(); ++k)
					for (int beta=0; beta < p_base->numberStringsOfType(k); ++beta) {
						if ((k==j) && (alpha==beta)) continue;

						for (int b=0; b< p_chain->stringLength(k); ++b) {
							complex<REAL> root_k = root (k, beta, b);
							deviation[a] *= sin (root_j - root_k + I* p_chain->zeta()) / sin (root_j - root_k - I*p_chain->zeta());
						}
					}
					// product over the same string
					int s1 = 2*(length_j-a);
					int s2 = 2*(length_j-a-1);
					int s3 = 2*a;
					int s4 = 2*(a-1);
					deviation[a] *= 1.0*(s1?sf_half_n_anis[s1]:1)*(s2?sf_half_n_anis[s2]:1)
						/( (s3?sf_half_n_anis[s3]:1)*(s4?sf_half_n_anis[s4]:1) );

					if (a>0) deviation[a] *= - deviation[a-1] ;
					sum_deviations += std::norm(deviation[a]);
				}
			}
		}
		return sqrt(sum_deviations);
	}


	/** longitudinal form factor |< GS | Sz | state >|^2 **/
	REAL Gap_State::longitudinalFormFactor (const State& state_mu) const
	{
		// different subspaces do not overlap
		if (state_mu.p_base->numberDown() != p_base->numberDown()) return 0.0;
		return 0.25*p_chain->length() * exp(Gap_State::lndetHm2P(state_mu));
	}


	/** transverse form factor |< GS | S- | state >|^2 **/
	REAL Gap_State::transverseFormFactor (const State& state_mu) const
	{
		// different subspaces do not overlap
		if (state_mu.p_base->numberDown() != p_base->numberDown()+1) return 0.0;
		return p_chain->length() * exp(Gap_State::lndetHminus(state_mu));
	}


	/** H minus 2 P determinant expression **/
	REAL Gap_State::lndetHm2P (const State& state_mu) const
	{
		const char* here = "Gap_State::lndetHm2P";
		if (state_mu.p_base->numberRoots() != p_base->numberRoots()) throw Exception (here, exc_NumberDown);

		for (int j=0; j< state_mu.p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < state_mu.p_base->numberStringsOfType(j); ++alpha) {
			if (( abs(state_mu.rapidity(j,alpha)+0.5*PI) < COLLISION_THRESHOLD ))
				throw Exception(here, exc_RapidityMinHalfPi, "in left state");
		}

		// check for collisions in right state lambda
		for (int j=0; j< p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha) {
			if (( abs(rapidity(j,alpha)+0.5*PI) < COLLISION_THRESHOLD ))
				throw Exception(here, exc_RapidityMinHalfPi, "in right state");
		}

		// collision checks
		/// collisions aren't supposed to happen any more, they should be caught by admissible()
		/// but aren't always...

		// check for collisions in left state mu
		for (int j=0; j< state_mu.p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < state_mu.p_base->numberStringsOfType(j); ++alpha) {
			for (int k=0; k< state_mu.p_base->numberTypes(); ++k)
			for (int beta=0; beta < state_mu.p_base->numberStringsOfType(k); ++beta) {
				if ((j==k) && (alpha==beta)) continue;
				if ( abs(rounderr((state_mu.rapidity(j,alpha) - state_mu.rapidity(k,beta))/PI)) < COLLISION_THRESHOLD )
					throw Exception(here, exc_RapidityCollision, "in left state");
			}
		}

		// check for collisions in right state lambda
		for (int j=0; j< p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha) {
			for (int k=0; k< p_base->numberTypes(); ++k)
			for (int beta=0; beta < p_base->numberStringsOfType(k); ++beta) {
				if ((j==k) && (alpha==beta)) continue;
				if ( abs(rounderr((rapidity(j,alpha) - rapidity(k,beta))/PI)) < COLLISION_THRESHOLD )
					throw Exception(here, exc_RapidityCollision, "in right state");
			}
		}


		// lambda, the right state of the matrix element
		int number_roots = p_base->numberRoots();
		complex<REAL>* root_lambda = new complex<REAL> [number_roots];
		complex<REAL>* sin_lambda = new complex<REAL> [number_roots];
		complex<REAL>* cos_lambda = new complex<REAL> [number_roots];
		for (int index=0, j=0; j< p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha)
		for (int a=0; a< p_chain->stringLength(j); ++a, ++index) {
			root_lambda[index] = root(j, alpha,a);
			sin_lambda[index] = sin(root_lambda[index]);
			cos_lambda[index] = cos(root_lambda[index]);
		}

		// mu, the left state of the matrix element
		complex<REAL>* root_mu = new complex<REAL> [number_roots];
		complex<REAL>* sin_mu = new complex<REAL> [number_roots];
		complex<REAL>* cos_mu = new complex<REAL> [number_roots];
		for (int index=0, j=0; j< state_mu.p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < state_mu.p_base->numberStringsOfType(j); ++alpha)
		for (int a=0; a< p_chain->stringLength(j); ++a, ++index) {
			root_mu[index] = state_mu.root(j, alpha, a);
			sin_mu[index] = sin(root_mu[index]);
			cos_mu[index] = cos(root_mu[index]);
		}

		// we only need these norm values at the very end, but calculate them now
		// so that there's no risk of having both the gaudin matrix
		// and hminus in memory at the same time.
		REAL norm_mu = state_mu.norm();
		REAL norm_lambda = norm();

		// these are all over the place.
		REAL cosh_zeta = cosh(p_chain->zeta());
		REAL sinh_zeta = sinh(p_chain->zeta());

		// calculate the F and G products
		complex<REAL>* ln_f0 = new complex<REAL> [number_roots];
		complex<REAL>* ln_f1 = new complex<REAL> [number_roots];
		complex<REAL>* ln_f2 = new complex<REAL> [number_roots];
		complex<REAL>* ln_g0 = new complex<REAL> [number_roots];
		complex<REAL>* ln_g1 = new complex<REAL> [number_roots];
		complex<REAL>* ln_g2 = new complex<REAL> [number_roots];
		REAL ln_pre = 0.0;
		for (int b=0, type_b=0; type_b < p_base->numberTypes(); ++type_b)
		for (int beta=0; beta < p_base->numberStringsOfType(type_b); ++beta) {
			int csi = b;
			int cs_length = p_chain->stringLength(type_b);
			for (int vertical_b = 0; vertical_b < cs_length; ++vertical_b, ++b) {
				ln_f0[b] = ln_f1[b] = ln_f2[b] = 0.0;
				ln_g0[b] = ln_g1[b] = ln_g2[b] = 0.0;

				complex<REAL> product_f0 = 1.0, product_f1 = 1.0, product_f2 = 1.0;
				complex<REAL> product_g0 = 1.0, product_g1 = 1.0, product_g2 = 1.0;
				complex<REAL> prefactor = 1.0;
				for (int k = 0; k < number_roots; ++k) {
					complex<REAL> s_lk_m_lb = sin_lambda[k] * cos_lambda[b] - cos_lambda[k] * sin_lambda[b];
					complex<REAL> c_lk_m_lb = cos_lambda[k] * cos_lambda[b] + sin_lambda[k] * sin_lambda[b];
					product_f0 *= s_lk_m_lb * cosh_zeta - I* c_lk_m_lb * sinh_zeta;
					if (k != b) product_f1 *= s_lk_m_lb;
					if ((k != b + 1) || (b-csi == cs_length-1))
						product_f2 *= s_lk_m_lb * cosh_zeta + I*c_lk_m_lb * sinh_zeta;

					complex<REAL> s_mk_m_lb = sin_mu[k] * cos_lambda[b] - cos_mu[k] * sin_lambda[b];
					complex<REAL> c_mk_m_lb = cos_mu[k] * cos_lambda[b] + sin_mu[k] * sin_lambda[b];
					product_g0 *= s_mk_m_lb * cosh_zeta - I* c_mk_m_lb * sinh_zeta;
					product_g1 *= s_mk_m_lb;
					product_g2 *= s_mk_m_lb * cosh_zeta + I* c_mk_m_lb * sinh_zeta;

					complex<REAL> s_ma_m_mk = sin_mu[b] * cos_mu[k] - cos_mu[b] * sin_mu[k];
					complex<REAL> c_ma_m_mk = cos_mu[b] * cos_mu[k] + sin_mu[b] * sin_mu[k];
					// filter out the zeroes in the denominator (need better implementation)
					if (::norm(s_ma_m_mk*cosh_zeta + I* c_ma_m_mk*sinh_zeta) > EQUAL_THRESHOLD)
						prefactor *= s_ma_m_mk*cosh_zeta + I* c_ma_m_mk*sinh_zeta;
					if (!(b+1==k) || (cs_length-1==vertical_b) )  {
						prefactor *= s_lk_m_lb*cosh_zeta + I*c_lk_m_lb *sinh_zeta;
					}

					// gather results in log (which avoids overflows)
					// but not too often as logs take time.
					if (!(k%50) || (k==number_roots-1)) {
						ln_f0[b] += log(product_f0);
						ln_f1[b] += log(product_f1);
						ln_f2[b] += log(product_f2);
						ln_g0[b] += log(product_g0);
						ln_g1[b] += log(product_g1);
						ln_g2[b] += log(product_g2);
						ln_pre -= log(abs( prefactor ));
						product_f0 = product_f1 = product_f2 = 1.0;
						product_g0 = product_g1 = product_g2 = 1.0;
						prefactor = 1.0;
					}
				}
			}
		}

		Square<complex<REAL> > matrix_hm2p (p_base->numberRoots());
		// a runs over all roots
		// type_a runs over all string types
		// alpha runs over all instances of a given string length
		// vertical_a runs over all roots in a given string instance.
		for (int a=0, type_a=0; type_a< state_mu.p_base->numberTypes(); ++type_a)
		for (int alpha=0; alpha < state_mu.p_base->numberStringsOfType(type_a); ++alpha)
		for (int vertical_a=0; vertical_a < state_mu.p_chain->stringLength(type_a); ++vertical_a, ++a) {

			for (int b=0, type_b=0; type_b< p_base->numberTypes(); ++type_b)
			for (int beta=0; beta < p_base->numberStringsOfType(type_b); ++beta) {
				int csi = b; // points to the start of the current string. inside root_lambda[]
				int cs_length = p_chain->stringLength(type_b);
				for (int vertical_b = 0; vertical_b < cs_length; ++vertical_b, ++b) {
					if ( (abs(imag(root_mu[a] - root_lambda[b])) < EQUAL_THRESHOLD)
						&& (abs(rounderr(real(root_mu[a] - root_lambda[b])/PI)) < EQUAL_THRESHOLD)
						) throw Exception (here, exc_Equal);
					if (1==cs_length) {
						// simple rapidities
						complex<REAL> s_ma_m_lcs = sin_mu[a] * cos_lambda[csi] - cos_mu[a] * sin_lambda[csi];
						complex<REAL> c_ma_m_lcs = cos_mu[a] * cos_lambda[csi] + sin_mu[a] * sin_lambda[csi];
						matrix_hm2p[a][b]=
							1.0/( ( s_ma_m_lcs * cosh_zeta - I* c_ma_m_lcs * sinh_zeta ) * s_ma_m_lcs )
							+ exp(ln_g2[csi] - ln_f2[csi] - ln_g0[csi] + ln_f0[csi])
							  / ( s_ma_m_lcs * (s_ma_m_lcs * cosh_zeta + I* c_ma_m_lcs * sinh_zeta))
							- 2.0 * exp(ln_f0[csi] - ln_g0[csi]) /(sq(sin_mu[a]) + sq(sf_half_n_anis[1]));
					}
					else if (vertical_b < cs_length-1) {
						// strings, element for "b<n"
						matrix_hm2p[a][b] = exp(-ln_g0[b])/( sin(root_mu[a] - root_lambda[b])
											* sin(root_mu[a] - root_lambda[b] + I*p_chain->zeta()) );
					}
					else {
						// strings, element for "b==n"
						// in these comments, subscripts indicate index as in paper; brackets indicate subscript as in this code.
						// lambda_1 == current_string_lambda[0]; I don't use lambda_0 and lambda_{n+1} since they're not real rapidities
						complex<REAL> sum_h (0.0, 0.0), sum_p (0.0, 0.0);

						// j=0 term ( a non-rapidity; remember, lambda_j with j==1 refers to the current root_lambda[b] )
						// paper: K * ( G_0 G_1 ) / ( F_0 F_1 )
						complex<REAL> s_ma_m_lcs = sin_mu[a] * cos_lambda[csi] - cos_mu[a] * sin_lambda[csi];
						complex<REAL> c_ma_m_lcs = cos_mu[a] * cos_lambda[csi] + sin_mu[a] * sin_lambda[csi];
						sum_h += exp(ln_g0[csi] + ln_g1[csi] - ln_f0[csi] - ln_f1[csi]
									- log( ( s_ma_m_lcs * cosh_zeta - I* c_ma_m_lcs * sinh_zeta ) * s_ma_m_lcs )	);

						// 0 < j < n terms: -d/d lambda K_j * ( G_j G_{j+1} ) / ( F_j F_{j+1} )
						for (int s_i = 0; s_i < cs_length-1; ++s_i) {
							complex<REAL> s_ma_m_lbpi = sin_mu[a] * cos_lambda[csi+s_i] - cos_mu[a] * sin_lambda[csi+s_i];
							complex<REAL> c_ma_m_lbpi = cos_mu[a] * cos_lambda[csi+s_i] + sin_mu[a] * sin_lambda[csi+s_i];
							sum_p += exp( ln_g1[csi+s_i] - ln_f1[csi+s_i] );
							sum_h -= exp( ln_g1[csi+s_i] + ln_g2[csi+s_i] - ln_f1[csi+s_i] - ln_f2[csi+s_i]
											+ log(	( c_ma_m_lbpi / s_ma_m_lbpi
													+ (c_ma_m_lbpi * cosh_zeta - I* s_ma_m_lbpi * sinh_zeta)
													/ (s_ma_m_lbpi * cosh_zeta + I* c_ma_m_lbpi * sinh_zeta)  )
												/	( s_ma_m_lbpi * (s_ma_m_lbpi * cosh_zeta + I* c_ma_m_lbpi * sinh_zeta)  )  )  );
						}

						// j=n term : H += K * ( G_n G_{n+1} ) / ( F_n F_{n+1} )
						complex<REAL> s_ma_m_lbpi = sin_mu[a] * cos_lambda[csi+cs_length-1] - cos_mu[a] * sin_lambda[csi+cs_length-1];
						complex<REAL> c_ma_m_lbpi = cos_mu[a] * cos_lambda[csi+cs_length-1] + sin_mu[a] * sin_lambda[csi+cs_length-1];
 						sum_p += exp( ln_g1[csi+cs_length-1] - ln_f1[csi+cs_length-1] );
						sum_h += exp( ln_g1[csi+cs_length-1] + ln_g2[csi+cs_length-1] - ln_f1[csi+cs_length-1] - ln_f2[csi+cs_length-1]
						 				- log( s_ma_m_lbpi * (s_ma_m_lbpi * cosh_zeta + I* c_ma_m_lbpi * sinh_zeta)) );

						complex<REAL> ln_normalisation =  ln_f0[csi] + ln_f1[csi] - ln_g0[b];		// N *= F1/G1 (now it is F0*F1/G1)
						for (int s_i = 1; s_i < cs_length; ++s_i) {
							ln_normalisation += ln_g1[csi+s_i]; 	//  N times product G_j, 2<=j<=n
						}
						ln_normalisation -= ln_g1[csi+cs_length-1];

						sum_p /= sq(sin_mu[a]) + sq(sf_half_n_anis[1]);
						matrix_hm2p[a][b] = (sum_h - 2.0* sum_p) * exp(ln_normalisation);
					} // end if

					if (!finite(matrix_hm2p[a][b])) {
						stringstream desc;
						desc << "on ("<<a<<", "<<b<<"); mu == "<<root_mu[a]<<"  lambda == "<<root_lambda[b];
						throw Exception(here, exc_NonFinite, desc.str());
					}
				}
			} // end lambda loop
		} // end mu loop

		REAL lndet_hm2p = lndetDestroy(matrix_hm2p);
		for (int b=0; b < number_roots; ++b)
			lndet_hm2p += real(ln_g0[b]);  // ln(abs(z)) == real(ln(z))

		REAL ln_form_factor = 0.0;
		for (int index=0; index < number_roots; ++index)
			ln_form_factor += 2.0 * log(abs(sin(root_mu[index] - 0.5*I*p_chain->zeta()))) - 2.0 * log(abs(sin(root_lambda[index] - 0.5*I*p_chain->zeta())));
		ln_form_factor += 2.0 * number_roots * log( sf_half_n_anis[2] );
		ln_form_factor += 2.0 * lndet_hm2p + ln_pre - norm_mu - norm_lambda;

		delete[] root_lambda; delete[] sin_lambda; delete[] cos_lambda;
		delete[] root_mu; delete[] sin_mu; delete[] cos_mu;
		delete[] ln_f0; delete[] ln_f1; delete[] ln_f2;
		delete[] ln_g0; delete[] ln_g1; delete[] ln_g2;
		return ln_form_factor;
	}


	/** H-minus determinant expression **/
	REAL Gap_State::lndetHminus (const State& state_mu) const
	{
		const char* here = "Gap_State::lndetHminus";
		if (state_mu.p_base->numberRoots() != p_base->numberRoots()+1) throw Exception (here, exc_NumberDown);

		for (int j=0; j< state_mu.p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < state_mu.p_base->numberStringsOfType(j); ++alpha) {
			if (( abs(state_mu.rapidity(j,alpha)+0.5*PI) < COLLISION_THRESHOLD ))
				throw Exception(here, exc_RapidityMinHalfPi, "in left state");
		}

		// check for collisions in right state lambda
		for (int j=0; j< p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha) {
			if (( abs(rapidity(j,alpha)+0.5*PI) < COLLISION_THRESHOLD ))
				throw Exception(here, exc_RapidityMinHalfPi, "in right state");
		}

		// collision checks
		/// collisions aren't supposed to happen any more, they should be caught by admissible()
		/// but aren't always...

		// check for collisions in left state mu
		for (int j=0; j< state_mu.p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < state_mu.p_base->numberStringsOfType(j); ++alpha) {
			for (int k=0; k< state_mu.p_base->numberTypes(); ++k)
			for (int beta=0; beta < state_mu.p_base->numberStringsOfType(k); ++beta) {
				if ((j==k) && (alpha==beta)) continue;
				if ( abs(rounderr((state_mu.rapidity(j,alpha) - state_mu.rapidity(k,beta))/PI)) < COLLISION_THRESHOLD )
					throw Exception(here, exc_RapidityCollision, "in left state");
			}
		}

		// check for collisions in right state lambda
		for (int j=0; j< p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha) {
			for (int k=0; k< p_base->numberTypes(); ++k)
			for (int beta=0; beta < p_base->numberStringsOfType(k); ++beta) {
				if ((j==k) && (alpha==beta)) continue;
				if ( abs(rounderr((rapidity(j,alpha) - rapidity(k,beta))/PI)) < COLLISION_THRESHOLD )
					throw Exception(here, exc_RapidityCollision, "in right state");
			}
		}

		// lambda, the right state of the matrix element
		int number_roots = p_base->numberRoots();
		complex<REAL>* root_lambda = new complex<REAL> [number_roots];
		complex<REAL>* sin_lambda = new complex<REAL> [number_roots];
		complex<REAL>* cos_lambda = new complex<REAL> [number_roots];
		for (int index=0, j=0; j< p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha)
		for (int a=0; a< p_chain->stringLength(j); ++a, ++index) {
			root_lambda[index] = root(j, alpha,a);
			sin_lambda[index] = sin(root_lambda[index]);
			cos_lambda[index] = cos(root_lambda[index]);
		}

		// mu, the left state of the matrix element
		int number_roots_mu = state_mu.p_base->numberRoots();
		complex<REAL>* root_mu = new complex<REAL> [number_roots_mu];
		complex<REAL>* sin_mu = new complex<REAL> [number_roots_mu];
		complex<REAL>* cos_mu = new complex<REAL> [number_roots_mu];
		for (int index=0, j=0; j< state_mu.p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < state_mu.p_base->numberStringsOfType(j); ++alpha)
		for (int a=0; a< p_chain->stringLength(j); ++a, ++index) {
			root_mu[index] = state_mu.root(j, alpha, a);
			sin_mu[index] = sin(root_mu[index]);
			cos_mu[index] = cos(root_mu[index]);
		}

		// we only need these norm values at the very end, but calculate them now
		// so that there's no risk of having both the gaudin matrix
		// and hminus in memory at the same time.
		REAL norm_mu = state_mu.norm();
		REAL norm_lambda = norm();

		// these are all over the place.
		REAL cosh_zeta = cosh(p_chain->zeta());
		REAL sinh_zeta = sinh(p_chain->zeta());

		// calculate the F and G products
		complex<REAL>* ln_f0 = new complex<REAL> [number_roots];
		complex<REAL>* ln_f1 = new complex<REAL> [number_roots];
		complex<REAL>* ln_f2 = new complex<REAL> [number_roots];
		complex<REAL>* ln_g0 = new complex<REAL> [number_roots];
		complex<REAL>* ln_g1 = new complex<REAL> [number_roots];
		complex<REAL>* ln_g2 = new complex<REAL> [number_roots];
		REAL ln_pre = 0.0;
		for (int b=0, type_b=0; type_b< p_base->numberTypes(); ++type_b)
		for (int beta=0; beta < p_base->numberStringsOfType(type_b); ++beta) {
			int csi = b;
			int cs_length = p_chain->stringLength(type_b);
			for (int vertical_b = 0; vertical_b < cs_length; ++vertical_b, ++b) {
				ln_f0[b] = ln_f1[b] = ln_f2[b] = 0.0;
				ln_g0[b] = ln_g1[b] = ln_g2[b] = 0.0;

				complex<REAL> product_f0 = 1.0, product_f1 = 1.0, product_f2 = 1.0;
				complex<REAL> product_g0 = 1.0, product_g1 = 1.0, product_g2 = 1.0;
				complex<REAL> prefactor = 1.0;
				for (int k = 0; k <= number_roots; ++k) {
					if (k < number_roots) {
						complex<REAL> s_lk_m_lb = sin_lambda[k] * cos_lambda[b] - cos_lambda[k] * sin_lambda[b];
						complex<REAL> c_lk_m_lb = cos_lambda[k] * cos_lambda[b] + sin_lambda[k] * sin_lambda[b];
						product_f0 *= s_lk_m_lb * cosh_zeta - I* c_lk_m_lb * sinh_zeta;
						if (k != b) product_f1 *= s_lk_m_lb;
						if ((k != b + 1) || (b-csi == cs_length-1))
							product_f2 *= s_lk_m_lb * cosh_zeta + I*c_lk_m_lb * sinh_zeta;

						if (!(b+1==k) || (cs_length-1==vertical_b) )  {
							prefactor *= s_lk_m_lb*cosh_zeta + I*c_lk_m_lb *sinh_zeta;
						}
					}

					complex<REAL> s_mk_m_lb = sin_mu[k] * cos_lambda[b] - cos_mu[k] * sin_lambda[b];
					complex<REAL> c_mk_m_lb = cos_mu[k] * cos_lambda[b] + sin_mu[k] * sin_lambda[b];
					product_g0 *= s_mk_m_lb * cosh_zeta - I* c_mk_m_lb * sinh_zeta;
					product_g1 *= s_mk_m_lb;
					product_g2 *= s_mk_m_lb * cosh_zeta + I* c_mk_m_lb * sinh_zeta;

					complex<REAL> s_ma_m_mk = sin_mu[b] * cos_mu[k] - cos_mu[b] * sin_mu[k];
					complex<REAL> c_ma_m_mk = cos_mu[b] * cos_mu[k] + sin_mu[b] * sin_mu[k];
					// filter out the zeroes in the denominator (need better implementation)
					if (::norm(s_ma_m_mk*cosh_zeta + I* c_ma_m_mk*sinh_zeta) > EQUAL_THRESHOLD)
						prefactor *= s_ma_m_mk*cosh_zeta + I* c_ma_m_mk*sinh_zeta;

					// gather results in log (which avoids overflows)
					// but not too often as logs take time.
					if (!(k%50) || (k==number_roots)) {
						ln_f0[b] += log(product_f0);
						ln_f1[b] += log(product_f1);
						ln_f2[b] += log(product_f2);
						ln_g0[b] += log(product_g0);
						ln_g1[b] += log(product_g1);
						ln_g2[b] += log(product_g2);
						ln_pre -= log(abs( prefactor ));
						product_f0 = product_f1 = product_f2 = 1.0;
						product_g0 = product_g1 = product_g2 = 1.0;
						prefactor = 1.0;
					}
				}
			}
		}

		// prefactor for last column b == M+1
		complex<REAL> prefactor = 1.0;
		for (int k = 0; k <= number_roots; ++k) {
			complex<REAL> s_ma_m_mk = sin_mu[number_roots_mu-1] * cos_mu[k] - cos_mu[number_roots_mu-1] * sin_mu[k];
			complex<REAL> c_ma_m_mk = cos_mu[number_roots_mu-1] * cos_mu[k] + sin_mu[number_roots_mu-1] * sin_mu[k];
			if (::norm(s_ma_m_mk*cosh_zeta + I* c_ma_m_mk*sinh_zeta) > EQUAL_THRESHOLD)
				prefactor *= s_ma_m_mk*cosh_zeta + I* c_ma_m_mk*sinh_zeta;
			// gather results in log
			if (!(k%50) || (k==number_roots)) {
				ln_pre -= log(abs( prefactor ));
				prefactor=1.0;
			}
		}

		// calculate the H- matrix
		Square< complex<REAL> > matrix_hminus (number_roots_mu);
		// type_a runs over all string types
		// alpha runs over all instances of a given string length
		// vertical_a runs over all roots in a given string instance.
		// a runs over all roots (first/row matrix index, goes with mu)
		for (int a=0, type_a=0; type_a< state_mu.p_base->numberTypes(); ++type_a)
		for (int alpha=0; alpha < state_mu.p_base->numberStringsOfType(type_a); ++alpha)
		for (int vertical_a=0; vertical_a < state_mu.p_chain->stringLength(type_a); ++vertical_a, ++a) {

			// last column of H-
			matrix_hminus[a][number_roots_mu-1] = 1.0 / (sq(sin_mu[a]) + sq(sf_half_n_anis[1]));

			// indices type_b, beta, vertical_B, b (second/column matrix index, goes with lambda): analogous to a.
			for (int b=0, type_b=0; type_b<p_base->numberTypes(); ++type_b)
			for (int beta=0; beta < p_base->numberStringsOfType(type_b); ++beta) {
				int csi = b;									// current string index: always points to the first root in the current string
				int cs_length = p_chain->stringLength(type_b); 	// length of the current string
				for (int vertical_b = 0; vertical_b < cs_length; ++vertical_b, ++b) {
					if ( (abs(imag(root_mu[a] - root_lambda[b])) < EQUAL_THRESHOLD)
						&& (abs(rounderr(real(root_mu[a] - root_lambda[b])/PI)) < EQUAL_THRESHOLD)
						) throw Exception (here, exc_Equal);
					if (1==cs_length) {
						// simple rapidities
						complex<REAL> s_ma_m_lcs = sin_mu[a] * cos_lambda[csi] - cos_mu[a] * sin_lambda[csi];
						complex<REAL> c_ma_m_lcs = cos_mu[a] * cos_lambda[csi] + sin_mu[a] * sin_lambda[csi];
						matrix_hminus[a][b] =
							1.0 / ( ( s_ma_m_lcs * cosh_zeta - I* c_ma_m_lcs * sinh_zeta ) * s_ma_m_lcs )
							+ exp(ln_g2[csi] - ln_f2[csi] - ln_g0[csi] + ln_f0[csi] )
							  / ( s_ma_m_lcs * (s_ma_m_lcs * cosh_zeta + I* c_ma_m_lcs * sinh_zeta))  ;
					}
					else if (vertical_b < cs_length-1)	{
						// strings, element for "b<n"
						matrix_hminus[a][b] = exp(-ln_g0[b])/( sin(root_mu[a] - root_lambda[b]) * sin(root_mu[a] - root_lambda[b] + I*p_chain->zeta()) );
					}
					else {
						// strings, element for "b==n"
						// in these comments, subscripts indicate index as in lateX notes; brackets indicate subscript as in this code.
						// lambda_1 == current_string_lambda[0] ; I don't use lambda_0 and lambda_{n+1} since they're not real rapidities
						complex<REAL> sum_h (0.0, 0.0);

						// j=0 term ( a non-rapidity; remember, lambda_j with j==1 refers to the current root_lambda[b] )
						// paper: K * ( G_0 G_1 ) / ( F_0 F_1 )
						complex<REAL> s_ma_m_lcs = sin_mu[a] * cos_lambda[csi] - cos_mu[a] * sin_lambda[csi];
						complex<REAL> c_ma_m_lcs = cos_mu[a] * cos_lambda[csi] + sin_mu[a] * sin_lambda[csi];
						sum_h += exp(ln_g0[csi] + ln_g1[csi] - ln_f0[csi] - ln_f1[csi]
								     - log( ( s_ma_m_lcs * cosh_zeta - I* c_ma_m_lcs * sinh_zeta ) * s_ma_m_lcs )	);

						// 0 < j < n terms: -d/d lambda K_j * ( G_j G_{j+1} ) / ( F_j F_{j+1} )
						for (int s_i = 0; s_i < cs_length-1; ++s_i) {
							complex<REAL> s_ma_m_lbpi = sin_mu[a] * cos_lambda[csi+s_i] - cos_mu[a] * sin_lambda[csi+s_i];
							complex<REAL> c_ma_m_lbpi = cos_mu[a] * cos_lambda[csi+s_i] + sin_mu[a] * sin_lambda[csi+s_i];
							sum_h -= exp(ln_g1[csi+s_i] + ln_g2[csi+s_i] - ln_f1[csi+s_i] - ln_f2[csi+s_i]
											+ log( (c_ma_m_lbpi / s_ma_m_lbpi
													+  (c_ma_m_lbpi * cosh_zeta - I* s_ma_m_lbpi * sinh_zeta)
													   / (s_ma_m_lbpi * cosh_zeta + I* c_ma_m_lbpi * sinh_zeta)	)
												/ ( s_ma_m_lbpi * (s_ma_m_lbpi * cosh_zeta + I* c_ma_m_lbpi * sinh_zeta)  )	)  );
						}

						// j=n term : K * ( G_n G_{n+1} ) / ( F_n F_{n+1} )
						complex<REAL> s_ma_m_lbpi = sin_mu[a] * cos_lambda[csi+cs_length-1] - cos_mu[a] * sin_lambda[csi+cs_length-1];
						complex<REAL> c_ma_m_lbpi = cos_mu[a] * cos_lambda[csi+cs_length-1] + sin_mu[a] * sin_lambda[csi+cs_length-1];
						sum_h += exp( ln_g1[csi+cs_length-1] + ln_g2[csi+cs_length-1] - ln_f1[csi+cs_length-1] - ln_f2[csi+cs_length-1]
										- log( s_ma_m_lbpi * (s_ma_m_lbpi * cosh_zeta + I* c_ma_m_lbpi * sinh_zeta)) );

						complex<REAL> ln_normalisation = ln_f0[csi] + ln_f1[csi] - ln_g0[b];
						for (int s_i = 1; s_i < cs_length; ++s_i)
							ln_normalisation += ln_g1[csi+s_i]; 	//  N times product G_j, 2<=j<=n
						ln_normalisation -= ln_g1[csi+cs_length-1];

						matrix_hminus[a][b] = sum_h * exp(ln_normalisation);
					} // end if
					if (!finite(matrix_hminus[a][b])) throw Exception(here, exc_NonFinite);
				}
			} // end lambda loop
		} // end mu loop

		REAL lndet_hminus = lndetDestroy (matrix_hminus);
		for (int b=0; b < number_roots; ++b)
			lndet_hminus += real(ln_g0[b]);  // ln(abs(z)) == real(ln(z))

		REAL ln_form_factor = 0.0;
		for (int index=0; index < p_base->numberRoots(); ++index)
			ln_form_factor += 2.0 * log(abs(sin(root_mu[index] - 0.5*I*p_chain->zeta()))) - 2.0 * log(abs(sin(root_lambda[index] - 0.5*I*p_chain->zeta())));
		ln_form_factor += 2.0 * log(abs(sin(root_mu[p_base->numberRoots()] - 0.5*I*p_chain->zeta())));
		ln_form_factor += (2.0 * number_roots + 2.0) * log( sf_half_n_anis[2] );
		ln_form_factor += 2.0 * lndet_hminus + ln_pre - norm_mu - norm_lambda;

		delete[] root_lambda; delete[] sin_lambda; delete[] cos_lambda;
		delete[] root_mu; delete[] sin_mu; delete[] cos_mu;
		delete[] ln_f0; delete[] ln_f1; delete[] ln_f2;
		delete[] ln_g0; delete[] ln_g1; delete[] ln_g2;
		return ln_form_factor;
	}


//}
