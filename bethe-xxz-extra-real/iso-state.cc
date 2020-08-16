#include "iso-state.h"

// used in matrices
#define EQUAL_THRESHOLD 1e-10
// to check if we have real part zero in quantum number
#define THRESHOLD_QUANTUM_NUMBER 1e-8

const char* exc_NotXXX = "cannot convert XXZ or gapped state to XXX";


	/** subspace ground state constructor **/
	XXX_State::XXX_State (const Base* const ground_base)
		: 	::State(ground_base)
	{
		//check consistency
		for (int i=1; i<p_base->numberTypes(); ++i)
			if (p_base->numberStringsOfType(i)) throw Exception ("XXX_State::XXX_State", exc_NotGroundBase);
		setFreeRapidities();
	}


	/** assignment **/
/*
	XXX_State& XXX_State::operator= (const XXX_State& rhs)
	{
		if (this != &rhs) State::operator= ( State(rhs) );
		return *this;
	}
*/



	/** roots of the bethe equation **/
	vector< complex<REAL> > XXX_State::roots (void) const
	{
		vector< complex<REAL> > root_vec (p_base->numberRoots());
		for (int index=0, j=0; j < p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha)
		for (int a=0; a < p_chain->stringLength(j); ++a, ++index) {
			root_vec[index] =  root(j, alpha, a); //rapidity(j, alpha) + 0.5 * ( p_chain->stringLength(j) - 1 - 2*a ) *I;
		}
		return root_vec;
	}


	/** energy, in units of J **/
	REAL XXX_State::calculateEnergy (void) const
	{
		its_energy = 0.0;
		for (int string_type=0; string_type< p_base->numberTypes(); ++string_type)
		for (int instance=0; instance < p_base->numberStringsOfType(string_type); ++instance)
		its_energy -= 2.0*(p_chain->stringLength(string_type))
				/ ( sq( 2.0* rapidity(string_type, instance) ) + sq(p_chain->stringLength(string_type)) );
		return its_energy;
	}




/** ***** ***** ******* ********** ************* ************ *********** *********** **/

	/** sum of scattering terms **/
	inline long double XXX_State::theta (const long double rap, const int length) const
	{	return 2.0 * atan( 2.0 * rap / (1.0*length) ) ; 	}

	/** bethe equation: right hand side **/
	long double XXX_State::rhs (const int j, const int alpha) const
	{
		long double scattering_term = 0.0;
		int length_j = p_chain->stringLength(j);

		for (int k=0; k < p_base->numberTypes(); ++k)
		for (int beta=0; beta < p_base->numberStringsOfType(k); ++beta) {
			int length_k = p_chain->stringLength(k);
			// epsilon of long double equals that of intermediate results
			long double lj_m_lk = rapidity(j, alpha) - rapidity(k, beta);
			scattering_term += XXX_State::theta (lj_m_lk, length_j + length_k);
			if (length_j != length_k)
				scattering_term += XXX_State::theta (lj_m_lk, abs(length_j-length_k));
			for (int i=1; i< min(length_j, length_k); ++i)
				scattering_term += 2.0 * XXX_State::theta (lj_m_lk, abs(length_j-length_k) + 2*i);
		}
		return PI*quantum_number(j, alpha) + scattering_term;
	}


	/** control-alt-delete. **/
	void XXX_State::setFreeRapidities (void)
	{
		// simply the bethe equation without scattering phases
		// j iterates over string types; alpha iterates over all instantiations of a string type; i.e. over all strings.
		for (int j=0; j < p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha)
			rapidity(j, alpha) = 0.5*p_chain->stringLength(j) * tan(0.5*PI*quantum_number(j, alpha)/(1.0*p_chain->length()));
		iterations=0;
		newton_iterations=0;
		convergence = NO_CONVERGENCE;
	}

	/** Bethe equation **/
	inline long double XXX_State::betheZero (const int j, const int alpha) const
	{
		return p_chain->length() * XXX_State::theta (rapidity(j, alpha), p_chain->string_length[j]) - XXX_State::rhs(j, alpha);
	};

	/** calculate Gaudin's matrix for this state and set the norm field **/

	/* d small theta/d lambda */
// 	inline long double XXX_State::thetaDerivative (const long double rap, const int length) const
// 	{	return 4.0 / ( 1.0*length  + 4.0*rap*rap/(1.0*length) ); 	}

	/* scattering derivative term (d big theta / d lambda) */
	long double XXX_State::scatteringDerivative (const int j, const int alpha, const int k, const int beta) const
	{
		int length_j = p_chain->stringLength(j);
		int length_k = p_chain->stringLength(k);
		//int prod_parity = p_chain->string_parity[j]*p_chain->string_parity[k];
		double lj_m_lk = rapidity(j, alpha) - rapidity(k, beta);
		long double big_theta = XXX_State::thetaDerivative (lj_m_lk, length_j+length_k);
		if (length_j != length_k)
			big_theta += XXX_State::thetaDerivative (lj_m_lk, abs(length_j - length_k));
		for (int i=1; i< min(length_j, length_k); ++i)
			big_theta += 2.0 * XXX_State::thetaDerivative (lj_m_lk, abs(length_j - length_k) + 2*i);

		return big_theta;
	}

	/* Gaudin's matrix */
	Square<REAL> XXX_State::matrixGaudin (void) const
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
					gaudin[index_j_alpha][index_k_beta] =
						p_chain->length() * XXX_State::thetaDerivative(rapidity(j, alpha), p_chain->stringLength(j));

					for (int l=0; l< p_base->numberTypes(); ++l)
					for (int gamma=0; gamma < p_base->numberStringsOfType(l); ++gamma) {
						if (rapidity.index(l, gamma) != index_j_alpha)
							gaudin[index_j_alpha][index_k_beta] -= XXX_State::scatteringDerivative (j, alpha, l, gamma);
					}
				}
				else 	gaudin[index_j_alpha][index_k_beta] = XXX_State::scatteringDerivative (j, alpha, k, beta);
				if (!finite(gaudin[index_j_alpha][index_k_beta])) throw Exception (here, exc_NonFinite);
			}
		}

		// set the norm value now that we're at it
		its_lnnorm = lndet(gaudin);
		return gaudin;
	}

	/** one iteration of Bethe's equation **/
	void XXX_State::iterate (void)
	{
		string here = "XXX_State::iterate";

		Strip<REAL> new_rapidity (p_base);
		REAL square_diff = 0.0;
		REAL one_over_n = 1.0/(1.0*p_chain->length());
		long double new_rap = 0.0;

		for (int j=0; j < p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha) {
			REAL current_rapidity = rapidity(j, alpha);	// this should be finite & <1
			for (int run = 0; run < run_max; ++run) {
				long double last_rap = new_rap;
				new_rap = 0.5*p_chain->stringLength(j) * tan( 0.5*one_over_n*XXX_State::rhs (j, alpha) );
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
		for (int j=0; j < p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha) {
			rapidity(j, alpha) = new_rapidity(j, alpha);
		}
		convergence = square_diff;
		++iterations;
	}

	/** deviation from string hypothesis: sum of squares of norm of deltas **/
	REAL XXX_State::stringDeviation (void) const
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
					deviation[a] = pow(  (root_j + 0.5*I) / (root_j - 0.5*I), - p_chain->length());

					// product over other strings
					for (int k=0; k < p_base->numberTypes(); ++k)
					for (int beta=0; beta < p_base->numberStringsOfType(k); ++beta) {
						if ((k==j) && (alpha==beta)) continue;
						for (int b=0; b< p_chain->stringLength(k); ++b) {
							complex<REAL> root_k = root (k, beta, b);
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

	/** longitudinal form factor squared, |< state_mu | Sz | *this >|^2 **/
	REAL XXX_State::longitudinalFormFactor (const State& state_mu) const
	{
		/// FIXME: this is only true if mu is the ground state.
		if (!admissible()) return 0.0;

		// different subspaces do not overlap
		if (state_mu.p_base->numberDown() != p_base->numberDown())
			return 0.0;
		if (p_base->numberInfiniteRapidities() == 0)
			return 0.25*p_chain->length() * exp(XXX_State::lndetHm2P(state_mu));
		// |<mu|Sz|lambda, infty>|^2 = |<mu|S-|lambda>|^2/(N-2(M-1))
		if (p_base->numberInfiniteRapidities() == 1)
			return exp(XXX_State::lndetHminus(state_mu)) / (1.0 - 2.0*(p_base->numberRoots())/(1.0*p_chain->length()));

		return 0.0;
	}


	/** transverse (spin lowering) form factor squared, |< state_mu | S- | (*this) >|^2 **/
	REAL XXX_State::transverseFormFactor (const State& state_mu) const
	{
		/// FIXME: this is only true if mu is the ground state.
		if (!admissible()) return 0.0;

		// different subspaces do not overlap. state_mu must have one more down spin since S- creates a down spin.
		if (state_mu.p_base->numberDown() != p_base->numberDown()+1)
			return 0.0;
		if (p_base->numberInfiniteRapidities())
			return 0.0; // |<mu|S-|lambda, infty>|^2 = 0

		return p_chain->length() * exp(XXX_State::lndetHminus(state_mu));
	}


	/** transverse (spin raising) form factor squared, |< state_mu | S+ | *this >|^2 **/
	REAL XXX_State::plusFormFactor (const State& state_mu) const
	{
		/// FIXME: this is only true if mu is the ground state.
		if (!admissible()) return 0.0;

		// different subspaces do not overlap
		int M_left = state_mu.p_base->numberDown();
		int M_right = p_base->numberDown();
		int N = p_chain->length();
		if (M_left != M_right-1)  return 0.0;

		if (p_base->numberInfiniteRapidities() == 0)
			return state_mu.transverseFormFactor (*this);
		// |<mu|S+|lambda, infty>|^2 = (4/N) |<mu|Sz|lambda>|^2 / (1-2M/N)
		// where M = M_{left_state} = M_{right_state} - 1 ; for #inf=1, this->p_base->number_roots == M_right_state-1.
		if (p_base->numberInfiniteRapidities() == 1)
			return exp(XXX_State::lndetHm2P(state_mu))/(1.0- 2.0*M_left/REAL(N));
		// where M-1 = M{left state}-1 ; numberRoots() is M{right_state}-2  and M{right_state} = M{left_state}+1
		// therefore M-1 == numberRoots()
		if (p_base->numberInfiniteRapidities() == 2)
  			return 	2.0 * N * exp(XXX_State::lndetHminus(state_mu)) / (  (1.0*N-2.0*M_left+2.0) * (1.0*N-2.0*M_left+1.0)  )  ;
  		// NOTE that the norm of <lambda|S+S+S-S-|lambda> is taken into account, so this is **NOT** (!) the composition of
  		// the S+|\lambda, \inf> -> Sz|\lambda>  and   Sz|\lambda,\inf> -> S- |\lambda>	formulas !!

  		// the above result is checked against exact diagonalisation.
		// simple composition would give the WRONG result
		//		return 	4.0 * N * exp(XXX_State::lndetHminus(state_mu)) / (  (1.0*N-2.0*M_left+2.0) * (1.0*N-2.0*M_left)  )  ;


		return 0.0;
	}


	/** H-2P determinant expression **/
	REAL  XXX_State::lndetHm2P (const State& state_mu) const
	{
		const char* here="XXX_State::lndetHm2P";

		// check for equal number of finite roots. this function doesn't know or care about infinite raps.
		if (state_mu.p_base->numberRoots() != p_base->numberRoots()) throw Exception (here, exc_NumberDown);

		// zero-sized matrix in lambda state: we have <|Sz(q, omega)|>^2 = N/4 det(H-2P) delta (q-q_mu) delta(omega-omega_mu) = N^2 delta (q-q_mu) delta(omega-omega_mu) / 4
		// |> is reference state |up,up,up,...>
		// NOTE that we need this for the S+- form factor
		//      but in case we'd like 'the connected Szz form factor for M=0' (which is zero) we would get the wrong result
		//      because this is an 'equal-rapidities' state
		// NOTE correction, we'd need all 'equal-rapidities' (connected) states for S+-, not only this one. for now, adapted sum rule instead.
		if (!p_base->numberRoots())
			throw Exception (here, exc_Equal, "(both states with no rapidities)");
			//return 0.0;


		// we only need these norm values at the very end, but calculate them now
		// so that there's no risk of having both the gaudin matrix
		// and hminus in memory at the same time.
		REAL norm_mu = state_mu.norm();
		REAL norm_lambda = norm();

		// lambda, the right state of the matrix element
		int number_roots = p_base->numberRoots();
		complex<REAL>* root_lambda = new complex<REAL> [number_roots];
		for (int index=0, j=0; j< p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha)
		for (int a=0; a< (j+1) ; ++a, ++index)
			root_lambda[index] = root(j, alpha, a);

		// mu, the left state of the matrix element
		int number_roots_mu = state_mu.p_base->numberRoots();
		complex<REAL>* root_mu = new complex<REAL> [number_roots_mu];
		for (int index=0, j=0; j< state_mu.p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < state_mu.p_base->numberStringsOfType(j); ++alpha)
		for (int a=0; a< (j+1) ; ++a, ++index)
			root_mu[index] = state_mu.root(j, alpha, a);

		// calculate the F and G products that are always needed, strings or no strings
		complex<REAL>* ln_f0 = new complex<REAL> [number_roots];
		complex<REAL>* ln_f1 = new complex<REAL> [number_roots];
		complex<REAL>* ln_f2 = new complex<REAL> [number_roots];
		complex<REAL>* ln_g0 = new complex<REAL> [number_roots];
		complex<REAL>* ln_g1 = new complex<REAL> [number_roots];
		complex<REAL>* ln_g2 = new complex<REAL> [number_roots];
		REAL ln_pre = 0.0;
		for (int b=0, type_b=0; type_b< p_base->numberTypes(); ++type_b)
		for (int beta=0; beta < p_base->numberStringsOfType(type_b); ++beta) {
			int csi = b;	// current string index: runs over the current string.
			int cs_length = (type_b+1);
			for (int vertical_b = 0; vertical_b < cs_length; ++vertical_b, ++b) {
				ln_f0[b] = ln_f1[b] = ln_f2[b] = 0.0;
				ln_g0[b] = ln_g1[b] = ln_g2[b] = 0.0;

				complex<REAL> product_f0 = 1.0, product_f1 = 1.0, product_f2 = 1.0;
				complex<REAL> product_g0 = 1.0, product_g1 = 1.0, product_g2 = 1.0;
				complex<REAL> prefactor = 1.0;
				for (int k = 0; k < number_roots; ++k) {
					product_f0 *= root_lambda[k] - root_lambda[b] - I;
					if (k != b) product_f1 *= root_lambda[k] - root_lambda[b];
					if ((k != b + 1) || (b-csi == cs_length-1))
						product_f2 *= root_lambda[k] - root_lambda[b] + I;

					product_g0 *= root_mu[k] - root_lambda[b] - I;
					product_g1 *= root_mu[k] - root_lambda[b];
					product_g2 *= root_mu[k] - root_lambda[b] + I;

					// these are the prefactors of the determinant expression,
					// but without the restriction that j != k, (which doesn't matter as we take absolute value)
					//
					// if we have strings, we have zero denominators  (root_lambda[b] - root_lambda[k] + I),   (i.e. root[b] == root[k] - I)
					// however, in that case, it cancels against a 1/0 in the unreduced Gaudin matrix.
					// since we use the reduced G. matrix, we should throw away these (near) zeroes


					if (!(b+1==k) || (cs_length-1==vertical_b) )  {
						prefactor *= root_lambda[k] - root_lambda[b] + I;
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
						ln_pre -= real(log( prefactor ));
						product_f0 = product_f1 = product_f2 = 1.0;
						product_g0 = product_g1 = product_g2 = 1.0;
						prefactor = 1.0;
					}
				}

				for (int c=0, type_c=0; type_c< state_mu.p_base->numberTypes(); ++type_c)
				for (int gamma=0; gamma < state_mu.p_base->numberStringsOfType(type_c); ++gamma) {
					int cs_length_mu = (type_c+1);
					for (int vertical_c = 0; vertical_c < cs_length_mu; ++vertical_c, ++c) {

						// these are the prefactors of the determinant expression,
						// but without the restriction that j != k, (which doesn't matter as we take absolute value)
						//
						// if we have strings, we have zero denominators  (root_mu[b] - root_mu[k] + I),   (i.e. root[b] == root[k] - I)
						// however, in that case, it cancels against a 1/0 in the unreduced Gaudin matrix.
						// since we use the reduced G. matrix, we should throw away these (near) zeroes
						//
						// this happens if 	b and k are subsequent roots in the same string, i.e.
						// *	b should be k+1
						// *    k and b should be in the same string, i.e. b is not the last in a string
						// was: ///if (::norm(root_mu[b] - root_mu[k] + I) > EQUAL_THRESHOLD)
						if ((c+1!=b) || (cs_length_mu-1==vertical_c) )  prefactor *= root_mu[b] - root_mu[c] + I;

						if (!(c%50) || (c==number_roots_mu-1)) {
							ln_pre -= real(log( prefactor ));
							prefactor = 1.0;
						}
					}
				}
			}
		}

		Square<complex<REAL> > matrix_hm2p (number_roots);
		// a runs over all roots
		// type_a runs over all string types
		// alpha runs over all instances of a given string length
		// vertical_a runs over all roots in a given string instance.
		for (int a=0, type_a=0; type_a< state_mu.p_base->numberTypes(); ++type_a)
		for (int alpha=0; alpha < state_mu.p_base->numberStringsOfType(type_a); ++alpha)
		for (int vertical_a=0; vertical_a < (type_a+1) ; ++vertical_a, ++a) {

			for (int b=0, type_b=0; type_b< p_base->numberTypes(); ++type_b)
			for (int beta=0; beta < p_base->numberStringsOfType(type_b); ++beta) {
				// points to the start of the current string. inside root_lambda[]
				int csi = b;
				int cs_length = (type_b+1);
				for (int vertical_b = 0; vertical_b < cs_length; ++vertical_b, ++b) {
					if (abs(root_mu[a] - root_lambda[b])<EQUAL_THRESHOLD ) throw Exception (here, exc_Equal);
					if (1==cs_length) {
						// simple rapidities
						matrix_hm2p[a][b] =
							// H
							1.0/( (root_mu[a] - root_lambda[csi] - I) * (root_mu[a] - root_lambda[csi]) )
							+ exp(ln_g2[csi] - ln_f2[csi] - ln_g0[csi] + ln_f0[csi])
								/ ( (root_mu[a] - root_lambda[csi] ) * (root_mu[a] - root_lambda[csi] + I)  )
							// minus 2P
							- 2.0* exp(ln_f0[csi] - ln_g0[csi]) / (sq(root_mu[a]) + 0.25);
					}
					else if (vertical_b < cs_length-1) 	{
						// strings, element for "b<n"
						matrix_hm2p[a][b] = exp(-ln_g0[b])/( (root_mu[a] - root_lambda[b]) * (root_mu[a] - root_lambda[b] + I) );
					}
					else {
						// strings, element for "b==n"
						// in these comments, subscripts indicate index as in paper; brackets indicate subscript as in this code.
						// lambda_1 == current_string_lambda[0]; I don't use lambda_0 and lambda_{n+1} since they're not real rapidities
						complex<REAL> sum_h (0.0, 0.0), sum_p (0.0, 0.0);

						// j=0 term: H += K * ( G_0 G_1 ) / ( F_0 F_1 )
						sum_h += exp(ln_g0[csi] + ln_g1[csi] - ln_f0[csi] - ln_f1[csi]
									- log( ( root_mu[a] - root_lambda[csi] - I) * (root_mu[a] - root_lambda[csi]) )	 );

						// 0 < j < n terms: H += -d/d lambda K_j * ( G_j G_{j+1} ) / ( F_j F_{j+1} )
						for (int s_i= 0; s_i< cs_length-1; ++s_i) {
							sum_p += exp( ln_g1[csi+s_i] - ln_f1[csi+s_i] );
							sum_h -= exp( ln_g1[csi+s_i] + ln_g2[csi+s_i] - ln_f1[csi+s_i]- ln_f2[csi+s_i]
										+ log(	( 1.0/(root_mu[a]- root_lambda[csi+s_i]) + 1.0/(root_mu[a]- root_lambda[csi+s_i] + I) )
											  / ( (root_mu[a] - root_lambda[csi+s_i])*(root_mu[a] - root_lambda[csi+s_i] + I) )  )    );
						}

						// j=n term : H += K * ( G_n G_{n+1} ) / ( F_n F_{n+1} )
						sum_p += exp( ln_g1[csi+cs_length-1] - ln_f1[csi+cs_length-1]);
						sum_h += exp( ln_g1[csi+cs_length-1] + ln_g2[csi+cs_length-1] - ln_f1[csi+cs_length-1] - ln_f2[csi+cs_length-1]
									- log( (root_mu[a] - root_lambda[csi+cs_length-1]) * (root_mu[a] - root_lambda[csi+cs_length-1]+ I)) );

						complex<REAL> ln_normalisation = ln_f0[csi] + ln_f1[csi] - ln_g0[b]; // equals N/G0
						for (int s_i = 1; s_i < cs_length; ++s_i)
							ln_normalisation += ln_g1[csi+s_i]; 				//  N times product G_j, 2<=j<=n
						ln_normalisation -= ln_g1[csi+cs_length-1]; //  N divided by G_n
						sum_p /= sq(root_mu[a]) + 0.25;
						matrix_hm2p[a][b] = (sum_h - 2.0* sum_p)*exp(ln_normalisation);
					} // end if

					if (!finite(matrix_hm2p[a][b])) {
						stringstream desc;
						desc << "on ("<<a<<", "<<b<<"); mu == "<<root_mu[a]<<"  lambda == "<<root_lambda[b];
						throw Exception(here, exc_NonFinite, desc.str());
					}
				}
			} // end lambda loop
		} // end mu loop

		REAL lndet_hm2p = lndetDestroy(matrix_hm2p); // actually ln(abs(det(..)))
		// apply the a-independent factors (i.e. we have divided the matrix columns by factors which we now apply to the determinant)
		for (int b=0; b < number_roots; ++b)
			lndet_hm2p += real(ln_g0[b]);  // ln(abs(z)) == real(ln(z))

		REAL ln_form_factor = 0.0;
		for (int index=0; index < p_base->numberRoots(); ++index) {
			ln_form_factor += 2.0 * log(abs(root_mu[index] - 0.5*I)) - 2.0 * log(abs(root_lambda[index] - 0.5*I));
		}
		ln_form_factor += 2.0 * lndet_hm2p + ln_pre - norm_mu - norm_lambda;

		delete[] root_lambda; delete[] root_mu;
		delete[] ln_f0; delete[] ln_f1; delete[] ln_f2;
		delete[] ln_g0; delete[] ln_g1; delete[] ln_g2;
		return ln_form_factor;

	}

	/** H minus determinant expression **/
	REAL XXX_State::lndetHminus (const State& state_mu) const
	{
		const char* here = "XXX_State::lndetHminus";

		if (state_mu.p_base->numberRoots() != p_base->numberRoots()+1) throw Exception (here, exc_NumberDown);

		// zero-sized matrix in lambda state: we have <mu|S-(q, omega)|>^2 = N det(H-) delta (q-q_mu) delta(omega-omega_mu) = delta (q-q_mu) delta(omega-omega_mu)
		// |> is reference state |up,up,up,...>
		if (!p_base->numberRoots()) return log(1.0/p_chain->length());

		// we only need these norm values at the very end, but calculate them now
		// so that there's no risk of having both the gaudin matrix
		// and hminus in memory at the same time.
		REAL norm_mu = state_mu.norm();
		REAL norm_lambda = norm();

		// lambda, the right state of the matrix element
		int number_roots = p_base->numberRoots();
		complex<REAL>* root_lambda = new complex<REAL> [number_roots];
		for (int index=0, j=0; j< p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha)
		for (int a=0; a< (j+1) ; ++a, ++index)
			root_lambda[index] = root(j, alpha, a);

		// mu, the left state of the matrix element
		int number_roots_mu = state_mu.p_base->numberRoots();
		complex<REAL>* root_mu = new complex<REAL> [number_roots_mu];
		for (int index=0, j=0; j< state_mu.p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < state_mu.p_base->numberStringsOfType(j); ++alpha)
		for (int a=0; a< (j+1) ; ++a, ++index)
			root_mu[index] = state_mu.root(j, alpha, a);

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
			int cs_length = (type_b+1);
			for (int vertical_b = 0; vertical_b < cs_length; ++vertical_b, ++b) {
				ln_f0[b] = ln_f1[b] = ln_f2[b] = 0.0;
				ln_g0[b] = ln_g1[b] = ln_g2[b] = 0.0;

				complex<REAL> product_f0 = 1.0, product_f1 = 1.0, product_f2 = 1.0;
				complex<REAL> product_g0 = 1.0, product_g1 = 1.0, product_g2 = 1.0;
				complex<REAL> prefactor = 1.0;
				for (int k = 0; k <= number_roots; ++k) {
					if (k < number_roots) {
						product_f0 *= root_lambda[k] - root_lambda[b] - I;
						if (k != b) product_f1 *= root_lambda[k] - root_lambda[b];
						if ((k != b + 1) || (b-csi == cs_length-1))
							product_f2 *= root_lambda[k] - root_lambda[b] + I;	// for 1-strings, always F2 == conj(F0)
						if (!(b+1==k) || (cs_length-1==vertical_b) )  {
							prefactor *= root_lambda[k] - root_lambda[b] + I;
						}
					}
					product_g0 *= root_mu[k] - root_lambda[b] -I;
					product_g1 *= root_mu[k] - root_lambda[b];
					product_g2 *= root_mu[k] - root_lambda[b] +I; // I don't use conj to allow complex raps


					// gather results in log (which avoids overflows)
					// but not too often as logs take time.
					if (!(k%50) || (k==number_roots)) {
						ln_f0[b] += log(product_f0);
						ln_f1[b] += log(product_f1);
						ln_f2[b] += log(product_f2);
						ln_g0[b] += log(product_g0);
						ln_g1[b] += log(product_g1);
						ln_g2[b] += log(product_g2);
						ln_pre -= real(log( prefactor ));			//TODO: is log(real(abs(.. faster?
						product_f0 = product_f1 = product_f2 = 1.0;
						product_g0 = product_g1 = product_g2 = 1.0;
						prefactor = 1.0;
					}
				}


				for (int c=0, type_c=0; type_c< state_mu.p_base->numberTypes(); ++type_c)
				for (int gamma=0; gamma < state_mu.p_base->numberStringsOfType(type_c); ++gamma) {
					int cs_length_mu = (type_c+1);
					for (int vertical_c = 0; vertical_c < cs_length_mu; ++vertical_c, ++c) {

						// these are the prefactors of the determinant expression,
						// but without the restriction that j != k, (which doesn't matter as we take absolute value)
						//
						// if we have strings, we have zero denominators  (root_mu[b] - root_mu[k] + I),   (i.e. root[b] == root[k] - I)
						// however, in that case, it cancels against a 1/0 in the unreduced Gaudin matrix.
						// since we use the reduced G. matrix, we should throw away these (near) zeroes
						//
						// this happens if 	b and k are subsequent roots in the same string, i.e.
						// *	b should be k+1
						// *    k and b should be in the same string, i.e. b is not the last in a string
						// was: ///if (::norm(root_mu[b] - root_mu[k] + I) > EQUAL_THRESHOLD)
						if ((c+1!=b) || (cs_length_mu-1==vertical_c) )  prefactor *= root_mu[b] - root_mu[c] + I;

						if (!(c%50) || (c==number_roots_mu-1)) {
							ln_pre -= real(log( prefactor ));
							prefactor = 1.0;
						}
					}
				}


			}
		}

		// prefactor for last column b == M+1
		complex<REAL> prefactor = 1.0;
		for (int c=0, type_c=0; type_c< state_mu.p_base->numberTypes(); ++type_c)
		for (int gamma=0; gamma < state_mu.p_base->numberStringsOfType(type_c); ++gamma) {
			int cs_length_mu = (type_c+1);
			for (int vertical_c = 0; vertical_c < cs_length_mu; ++vertical_c, ++c) {

				// these are the prefactors of the determinant expression,
				// but without the restriction that j != k, (which doesn't matter as we take absolute value)
				//
				// if we have strings, we have zero denominators  (root_mu[b] - root_mu[k] + I),   (i.e. root[b] == root[k] - I)
				// however, in that case, it cancels against a 1/0 in the unreduced Gaudin matrix.
				// since we use the reduced G. matrix, we should throw away these (near) zeroes
				//
				// this happens if 	b and k are subsequent roots in the same string, i.e.
				// *	b should be k+1
				// *    k and b should be in the same string, i.e. b is not the last in a string
				// was: ///if (::norm(root_mu[b] - root_mu[k] + I) > EQUAL_THRESHOLD)
				if ((c+1!=number_roots_mu-1) || (cs_length_mu-1==vertical_c) )  prefactor *= root_mu[number_roots_mu-1] - root_mu[c] + I;

				if (!(c%50) || (c==number_roots_mu-1)) {
					ln_pre -= real(log( prefactor ));
					prefactor = 1.0;
				}
			}
		}

		Square< complex<REAL> > matrix_hminus (number_roots_mu);
		// a runs over all roots
		// type_a runs over all string types
		// alpha runs over all instances of a given string length
		// vertical_a runs over all roots in a given string instance.
		for (int a=0, type_a=0; type_a< state_mu.p_base->numberTypes(); ++type_a)
		for (int alpha=0; alpha < state_mu.p_base->numberStringsOfType(type_a); ++alpha)
		for (int vertical_a=0; vertical_a < (type_a+1) ; ++vertical_a, ++a) {

			// last column of H-
			matrix_hminus[a][number_roots] =	1.0 / (sq(root_mu[a]) + 0.25);

			for (int b=0, type_b=0; type_b<p_base->numberTypes(); ++type_b)
			for (int beta=0; beta < p_base->numberStringsOfType(type_b); ++beta) {
				int csi = b;				  // points to the start of the current string. inside root_lambda[]
				int cs_length =  (type_b+1) ; // current string length
				for (int vertical_b = 0; vertical_b < cs_length; ++vertical_b, ++b) {
					if (abs(root_mu[a] - root_lambda[b]) < EQUAL_THRESHOLD) throw Exception (here, exc_Equal);
					if (1==cs_length) {
						// simple rapidities
						matrix_hminus[a][b] =
							1.0/( (root_mu[a] - root_lambda[csi] - I) * (root_mu[a] - root_lambda[csi]) )
							+ exp( ln_g2[csi] - ln_f2[csi]- ln_g0[csi] + ln_f0[csi] )
							  / ( (root_mu[a] - root_lambda[csi] ) * (root_mu[a] - root_lambda[csi] + I)  );
					}
					else if (vertical_b < cs_length-1) {
						// strings, element for "b<n"
						matrix_hminus[a][b] = exp(-ln_g0[b]) /( (root_mu[a] - root_lambda[b]) * (root_mu[a] - root_lambda[b] + I) );
					}
					else {
						// strings, element for "b==n"
						// in these comments, subscripts indicate index as in paper; brackets indicate subscript as in this code.
						// lambda_1 == current_string_lambda[0]; I don't use lambda_0 and lambda_{n+1} since they're not real rapidities
						complex<REAL> sum_h (0.0, 0.0);

						// j=0 term ( a non-rapidity; remember, lambda_j with j==1 refers to the current root_lambda[b] )
						// paper: K * ( G_0 G_1 ) / ( F_0 F_1 )
						sum_h += 	exp( ln_g0[csi] + ln_g1[csi] - ln_f0[csi] - ln_f1[csi]
										- log( (root_mu[a] - root_lambda[csi] - I) * (root_mu[a] - root_lambda[csi]) ) );

						// 0 < j < n terms: -d/d lambda K_j * ( G_j G_{j+1} ) / ( F_j F_{j+1} )
						for (int s_i=0; s_i < cs_length-1 ; ++s_i) { // s_i: string index
							sum_h -= exp( ln_g1[csi+s_i] + ln_g2[csi+s_i] - ln_f1[csi+s_i] - ln_f2[csi+s_i]
										+ log(	( 1.0/(root_mu[a] - root_lambda[csi+s_i]) + 1.0/(root_mu[a] - root_lambda[csi+s_i] +I) )
												/	( (root_mu[a] - root_lambda[csi+s_i]) * (root_mu[a] - root_lambda[csi+s_i] +I) )  )	);
						}

						// j=n term : K * ( G_n G_{n+1} ) / ( F_n F_{n+1} )
						sum_h += exp( ln_g1[csi+cs_length-1] + ln_g2[csi+cs_length-1] - ln_f1[csi+cs_length-1] - ln_f2[csi+cs_length-1]
									- log( (root_mu[a] - root_lambda[csi+cs_length-1]) * (root_mu[a] - root_lambda[csi+cs_length-1]+ I)) );

						complex<REAL> ln_normalisation = ln_f1[csi] + ln_f0[csi] - ln_g0[b]; // equals N /G0
						for (int s_i = 1; s_i < cs_length; ++s_i)
							ln_normalisation += ln_g1[csi+s_i]; 				//  N times product G_j, 2<=j<=n
						ln_normalisation -= ln_g1[csi+cs_length-1]; //  N divided by G_n
						matrix_hminus[a][b] = sum_h * exp(ln_normalisation);
					} // end if

					if (!finite(matrix_hminus[a][b])) {
						stringstream desc;
						desc << "on ("<<a<<", "<<b<<"); mu == "<<root_mu[a]<<"  lambda == "<<root_lambda[b];
						throw Exception(here, exc_NonFinite, desc.str());
					}
				}
			} // end lambda loop
		} // end mu loop

		REAL lndet_hminus = lndetDestroy (matrix_hminus); // actually ln(abs(det(..))). matrix_hminus gets destroyed.
		// apply the a-independent factors (i.e. we have divided the matrix columns by factors which we now apply to the determinant)
		for (int b=0; b < number_roots; ++b)
			lndet_hminus += real(ln_g0[b]);  // ln(abs(z)) == real(ln(z))

		REAL ln_form_factor=0.0;
		for (int index=0; index < p_base->numberRoots(); ++index)
			ln_form_factor += 2.0 * log(abs( (root_mu[index] - 0.5*I) )) - 2.0 * log(abs( (root_lambda[index] - 0.5*I) ));
		ln_form_factor += 2.0 * log(abs( (root_mu[state_mu.p_base->numberRoots()-1] - 0.5*I) ));
		ln_form_factor += 2.0 * lndet_hminus + ln_pre - norm_mu - norm_lambda;

		delete[] ln_f0; delete[] ln_f1; delete[] ln_f2;
		delete[] ln_g0; delete[] ln_g1; delete[] ln_g2;
		delete[] root_lambda; delete[] root_mu;
		return ln_form_factor;
	}


	/** calculate the quantum numbers (2I) as they are in the complex (non-strung) Bethe equations, from I$ = I1+I2+N/2, I1-I2 = [wide] **/
	vector< int > XXX_State::calculateBetheQuantumNumbers (void) const
	{
		const int twostring_type=1;

		// currently, we ignore infinite rapidities
		vector< int > bethe_quantum_number (p_base->numberRoots());
		for (int index_j_alpha = 0, j=0; j< p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha)
		for (int a=0; a< (j+1) ; ++a, ++index_j_alpha) {
			if (0==j) bethe_quantum_number[index_j_alpha] = quantum_number(j, alpha);
			else if (twostring_type==j) {
				//bool narrow = ( p_chain->length()/2 - p_base->numberRoots() + quantum_number[twostring_type][0]/2 )%2;
					// criterion for narrow/wide pair, from I$ == I1 + I2 + N/2 ; I1 - I2 = [wide]
				bool narrow = ( (-p_chain->length() - p_base->numberStringsOfType(twostring_type)+1 + quantum_number(twostring_type, alpha))/2 + p_base->numberRoots() )%2;
				bethe_quantum_number[index_j_alpha] = quantum_number(j, alpha)/2 + p_chain->length()/2 + (narrow?0:(a?1:-1));
			}
			else throw Exception ("XXX_State::calculateBetheQuantumNumbers", "not implemented for strings longer than 2");
		}
		return bethe_quantum_number;
	}




	/** calculate the quantum numbers (J) as they are in the complex (non-strung) Bethe equations, from the solution **/
	vector< complex<REAL> > XXX_State::calculateBetheI (void) const
	{
		// currently, we ignore infinite rapidities
		vector< complex<REAL> > bethe_roots = roots();
		vector< complex<REAL> > bethe_i (p_base->numberRoots());
		for (int alpha=0; alpha < p_base->numberRoots(); ++alpha) {
			if (abs(real(bethe_roots[alpha])) < THRESHOLD_QUANTUM_NUMBER && abs(abs(imag(bethe_roots[alpha]))-0.5)< THRESHOLD_QUANTUM_NUMBER ) {
				// numberRoots is number of finite roots. if even, there is no zero; if odd, there is a zero.
				int root_sign = isgn(imag(bethe_roots[alpha]));
				if (p_base->numberRoots()%2)
					// odd. there's a zero q.n.
					bethe_i[alpha] = 0.25*(p_chain->length() - root_sign*(p_chain->length()%4));
				else
					// even. no zero.
					bethe_i[alpha] = 0.25*(p_chain->length() - root_sign*(2 - p_chain->length()%4));
				continue;
			}
			complex<REAL> lhs = atan(2.0*bethe_roots[alpha]) * REAL(p_chain->length());
			for (int beta=0; beta<p_base->numberRoots(); ++beta){
				if (alpha!=beta) lhs -= atan(bethe_roots[alpha] - bethe_roots[beta]);
			}
			bethe_i[alpha] = lhs/PI;
		}
		return bethe_i;
	}

	/** calculate the quantum numbers (J) as they are in the complex (non-strung) Bethe equations, from the solution. round to int **/
	///FIXME: the naming of these three methods is horrible
	vector< int > XXX_State::calculateBethe2I (void) const
	{
		vector< complex<REAL> > bethe_i = calculateBetheI();
		vector<int> bethe_quantum_number (bethe_i.size());
		for (int i=0; i<bethe_i.size();++i)
			bethe_quantum_number[i] = int( round(2.0*real(bethe_i[i])) );
		return bethe_quantum_number;
	}

	/** check the complex Bethe equations **/
	REAL XXX_State::betheError(void) const
	{
		vector<complex<REAL> > bethe_j = calculateBetheI();
		REAL error = 0.0;
		for (int i=0; i<bethe_j.size();++i)
			error += ::norm(bethe_j[i] -  0.5*round(2.0*real(bethe_j[i])));
		return error;
	/*
		REAL error = 0.0;
		vector< complex<REAL> > bethe_roots = roots();
		for (int alpha=0; alpha < p_base->numberRoots(); ++alpha) {
			if (	abs(real(bethe_roots[alpha])) < THRESHOLD_QUANTUM_NUMBER
					&& abs(abs(imag(bethe_roots[alpha]))-0.5)< THRESHOLD_QUANTUM_NUMBER )
				continue;
			complex<REAL> lhs = pow( (bethe_roots[alpha]-0.5*I)/(bethe_roots[alpha]+0.5*I), p_chain->length());
			for (int beta=0; beta<p_base->numberRoots(); ++beta){
				if (alpha!=beta) lhs /= (bethe_roots[alpha]-bethe_roots[beta]-I)/(bethe_roots[alpha]-bethe_roots[beta]+I);
			}
			error += ::norm(1.0 - lhs);
		}
		return error;
	*/
	}


