#include "xxz-state.h"

#define EQUAL_THRESHOLD 1e-10
const char* exc_NotXXZ = "cannot convert XXX or gapped state to XXZ";



	/** subspace ground state constructor **/
	XXZ_State::XXZ_State (const Base* const ground_base)
		: 	::State(ground_base), tfh_rapidity(ground_base)
	{
		//check consistency
		for (int i=1; i<p_base->numberTypes(); ++i){
			if (p_base->numberStringsOfType(i)) throw Exception ("XXZ_State::XXZ_State", exc_NotGroundBase);
		}
		setTables();
		setFreeRapidities();
	}




	/** assignment **/
	/*
	XXZ_State& XXZ_State::operator= (const XXZ_State& rhs)
	{
		if (this != &rhs) {
			State::operator= ( State(rhs) );
			tfh_rapidity = rhs.tfh_rapidity;
			setTables();
		}
		return *this;
	}
	*/


	/** copy from generic state **/
	XXZ_State::XXZ_State (const State& original): 	State(original), tfh_rapidity(original.p_base)
	{
		if (original.p_chain->delta() <= 0.0 || original.p_chain->delta() >= 1.0)
			throw Exception("XXZ_State::XXZ_State", exc_NotXXZ);
		// set tanh table
		for (int e=0; e<rapidity.numberElements(); ++e)
			setRapidityAtIndex(e, original.rapidity.element(e));
		setTables();
	};


	/** set tables **/
	void XXZ_State::setTables(void)
	{
		// last string is longest string
		int table_size = 2* p_chain->stringLength(p_chain->numberTypes() -1);
		tf_half_n_anis = new REAL[ table_size ];
		sf_half_n_anis = new REAL[ table_size ];
		for (int n=0; n < table_size; ++n) {
			tf_half_n_anis[n] = tan (0.5*n*p_chain->zeta());
			sf_half_n_anis[n] = sin (0.5*n*p_chain->zeta());
		}
	}

	/** set rapidities and related data fields **/
	/* from given rapidity */
	void XXZ_State::setRapidity (const int j, const int alpha, const long double new_rapidity)
	{
		rapidity(j, alpha) = new_rapidity;
		tfh_rapidity(j, alpha) = tanh(new_rapidity);
	}

	void XXZ_State::setRapidityAtIndex (const int index, const long double new_rapidity)
	{
		rapidity.element(index) = new_rapidity;
		tfh_rapidity.element(index) = tanh(new_rapidity);
	}

	/** roots of the bethe equation **/
	/* full vector */
	vector< complex<REAL> > XXZ_State::roots (void) const
	{
		vector< complex<REAL> > root (p_base->numberRoots());
		for (int index=0, j=0; j < p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha)
		for (int a=0; a < p_chain->stringLength(j); ++a, ++index) {
			root[index] =  rapidity(j, alpha)
				+ 0.5 * I * (p_chain->zeta() * (p_chain->string_length[j] - 1 - 2*a) )
				+ 0.25* I * (PI*(1-p_chain->string_parity[j]));
		}
		return root;
	}

	/* only one root */
	complex<REAL> XXZ_State::root (const int j, const int alpha, const int a) const
	{
		return	rapidity(j, alpha)
					+ 0.5 * I * (p_chain->zeta() * (p_chain->string_length[j] - 1 - 2*a) )
					+ 0.25* I * (PI*(1-p_chain->string_parity[j]));
	}


	/** energy, in units of J **/
	REAL XXZ_State::calculateEnergy (void) const
	{
		its_energy = 0.0;
		for (int string_type=0; string_type< p_base->numberTypes(); ++string_type)
			for (int instance=0; instance<p_base->numberStringsOfType(string_type); ++instance)
				its_energy -=	sin(p_chain->zeta()) * sin( p_chain->string_length[string_type] * p_chain->zeta())
								/
								(	p_chain->string_parity[string_type] * cosh( 2.0* rapidity(string_type, instance) )
									- cos(  p_chain->string_length[string_type] *p_chain->zeta() )
								);
		return its_energy;
	}


	/** internal methods, for readability **/
	/* d small theta/d lambda */
	inline long double XXZ_State::thetaDerivative (const long double rapidity, const int length, const int parity) const
 	{
		return 	(parity==1)	?sf_half_n_anis[2*length]/(sq(sinh(rapidity)) + sq(sf_half_n_anis[length]))
							:sf_half_n_anis[2*length]/(-sq(cosh(rapidity)) + sq(sf_half_n_anis[length]));
	}

	/* small theta  */
	inline long double XXZ_State::thetaTfh(const long double th_rapidity, const int length, const int parity) const
	{	return (parity==1)
				? 2.0 * atan( th_rapidity / tf_half_n_anis[length] )
				: -2.0 * atan( th_rapidity * tf_half_n_anis[length] ) ;
	}

	/* rapidity from phase */
	inline long double XXZ_State::tfhInvTheta (const long double phase, const int string_type) const
	{
		return  tan(-0.5*phase) / tan( 0.5*p_chain->stringLength(string_type)*p_chain->zeta() - 0.25*PI*(1.0+p_chain->string_parity[string_type]) ) ;
	}


	/** bethe equation: right hand side **/
	long double XXZ_State::rhs (const int j, const int alpha) const
	{
		long double scattering_term = 0.0;
		int length_j = p_chain->stringLength(j);
		int parity_j = p_chain->string_parity[j];
		REAL current_tanh_rapidity = tfh_rapidity(j, alpha);

		for (int k=0; k < p_base->numberTypes(); ++k)
		for (int beta=0; beta < p_base->numberStringsOfType(k); ++beta) {
			int length_k = p_chain->string_length[k];
			int prod_parity = parity_j*p_chain->string_parity[k];
			// epsilon of long double equals that of intermediate results
			long double th_lj_m_lk =
					(current_tanh_rapidity - tfh_rapidity(k, beta))
					/ (1.0 - current_tanh_rapidity*tfh_rapidity(k, beta));
			scattering_term += XXZ_State::thetaTfh (th_lj_m_lk, length_j + length_k, prod_parity);
			if (length_j != length_k)
				scattering_term += XXZ_State::thetaTfh (th_lj_m_lk, abs(length_j-length_k), prod_parity);
			for (int i=1; i< min(length_j, length_k); ++i)
				scattering_term += 2.0 * XXZ_State::thetaTfh (th_lj_m_lk, abs(length_j-length_k) + 2*i, prod_parity);
		}
		return PI*quantum_number(j, alpha) + scattering_term;
	}


	/** bethe equation: left hand side **/
	inline long double XXZ_State::betheZero (const int j, const int alpha) const
	{
		return
			p_chain->length() * XXZ_State::thetaTfh (tfh_rapidity(j, alpha), p_chain->stringLength(j), p_chain->string_parity[j])
				- XXZ_State::rhs(j, alpha);
	}



	/** initial state for iteration **/
	void XXZ_State::setFreeRapidities (void)
	{
		// simply the bethe equation without scattering phases
		// j iterates over string types; alpha iterates over all instantiations of a string type; i.e. over all strings.

		for (int j=0; j < p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha) {
			tfh_rapidity(j, alpha) = XXZ_State::tfhInvTheta ( PI*quantum_number(j, alpha)/(1.0*p_chain->length()), j);
			rapidity(j, alpha) = atanh(tfh_rapidity(j, alpha));
			if (!finite(rapidity(j, alpha))) {
				rapidity(j, alpha) = 0.0;
				tfh_rapidity(j, alpha) = 0.0;
			}
		}

		iterations=0;
		newton_iterations=0;
		convergence = NO_CONVERGENCE;

	}



/** ***** ***** ******* ********** ************* ************ *********** *********** **/


	/** Gaudin's matrix **/
	/* scattering derivative term (d big theta / d lambda) */
	long double XXZ_State::scatteringDerivative (const int j, const int alpha, const int k, const int beta) const
	{
		int length_j = p_chain->stringLength(j);
		int length_k = p_chain->stringLength(k);
		int prod_parity = p_chain->string_parity[j]*p_chain->string_parity[k];
		double lj_m_lk = rapidity(j, alpha) - rapidity(k, beta);
		long double big_theta = XXZ_State::thetaDerivative (lj_m_lk, length_j+length_k, prod_parity );
		if (length_j != length_k)
			big_theta += XXZ_State::thetaDerivative (lj_m_lk, abs(length_j - length_k), prod_parity);
		for (int i=1; i< min(length_j, length_k); ++i)
			big_theta += 2.0 * XXZ_State::thetaDerivative (lj_m_lk, abs(length_j - length_k) + 2*i, prod_parity);

		return big_theta;
	}

	/* derivative of left hand side of bethe equation */
	inline long double XXZ_State::lhsDerivative(const int j, const int alpha) const
	{
		return p_chain->length() * XXZ_State::thetaDerivative(rapidity(j, alpha), p_chain->stringLength(j), p_chain->string_parity[j]);
	}

	/* Gaudin's matrix */
	Square<REAL> XXZ_State::matrixGaudin (void) const
	{
		const char* here = "XXZ_State::matrixGaudin";
		Square<REAL> gaudin (rapidity.numberElements());
		// size is the number of string types in rapidity; numberElements is the number of strings in rapidity.

		for (int j=0; j< p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha) {

			int index_j_alpha = rapidity.index(j, alpha);
			for (int k=0; k< p_base->numberTypes(); ++k)
			for (int beta=0; beta < p_base->numberStringsOfType(k); ++beta) {

				int index_k_beta = rapidity.index(k, beta);
				if (index_j_alpha == index_k_beta) {
					gaudin[index_j_alpha][index_k_beta] = lhsDerivative (j, alpha);

					for (int l=0; l< p_base->numberTypes(); ++l)
					for (int gamma=0; gamma < p_base->numberStringsOfType(l); ++gamma) {
						if (rapidity.index(l, gamma) != index_j_alpha)
							gaudin[index_j_alpha][index_k_beta] -= scatteringDerivative (j, alpha, l, gamma);
					}
				}
				else 	gaudin[index_j_alpha][index_k_beta] = scatteringDerivative (j, alpha, k, beta);
				if (!finite(gaudin[index_j_alpha][index_k_beta])) throw Exception (here, exc_NonFinite);
			}
		}

		// set the norm value now that we're at it
		its_lnnorm = lndet(gaudin);
		return gaudin;
	}


	/** one iteration of Bethe's equation **/

	void XXZ_State::iterate (void)
	{
		string here = "State::iterate";
		Strip<REAL> new_rapidity (p_base);
		Strip<REAL> new_tfh_rapidity (p_base);
		REAL square_diff = 0.0;
		REAL one_over_n = 1.0/(1.0*p_chain->length());
		long double new_tfh_rap = 1.0;

		for (int j=0; j < p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha) {
			REAL current_tfh_rapidity = tfh_rapidity(j, alpha);	// this should be finite & <1
			REAL run_finish = sgn(current_tfh_rapidity)? 1.0*sgn(current_tfh_rapidity) :1.0  ;
			for (int run = 0; run < run_max; ++run) {
				new_tfh_rap = XXZ_State::tfhInvTheta (one_over_n * XXZ_State::rhs (j, alpha), j);

				if (finite(new_tfh_rap) && (abs(new_tfh_rap) < 1.0)) break;
				if (run == run_max/2) {
					// we're halfway & not converged: try the other way around.
					tfh_rapidity(j, alpha)= current_tfh_rapidity;
					run_finish = - run_finish;
				}
				// if it doesn't work, slowly raise the offending rapidity and see if we can get something reasonable.
				// NOTE: since rhs only uses tfh_rapidity, we don't calculate the atan.
				// thus rapidity[][] is temporarily out of sync !!
				tfh_rapidity(j, alpha) = run_finish * (1.0-run_sloth) + tfh_rapidity(j, alpha) * run_sloth;
				if (abs(tfh_rapidity(j, alpha)) == 1.0) break;
			}
			tfh_rapidity(j, alpha) = current_tfh_rapidity; // reset the old value, don't mess with current rapidities.
			new_tfh_rapidity(j, alpha) = new_tfh_rap;
			new_rapidity(j, alpha) = atanh(new_tfh_rapidity(j, alpha));
			if (!finite(new_rapidity(j, alpha)) || (abs(new_tfh_rap) >= 1.0)) throw Exception (here, exc_Runaway);
			square_diff += sq(rapidity(j, alpha) - new_rapidity(j, alpha));
		}
		// if we have NaN square diff, throw exception. *this shouldn't happen.*
		if  (!finite(square_diff)) throw Exception (here, exc_NonFinite);

		// set the new values for the rapidities
		for (int j=0; j < p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha) {
			rapidity(j, alpha) = new_rapidity(j, alpha);
			tfh_rapidity(j, alpha) = new_tfh_rapidity(j, alpha);
			// avoid false positive convergence: we can't trust the arctanh for number this large
			if (fabs(new_rapidity(j, alpha)) > ATANH_THRESHOLD) throw Exception (here, exc_Uncontrol);
		}
		convergence = square_diff;
		++iterations;
	}


	/** deviation from string hypothesis: sum of squares of norm of deltas **/
	REAL XXZ_State::stringDeviation (void) const
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
					deviation[a] = pow(  sinh(root_j + 0.5*I*p_chain->zeta()) / sinh(root_j - 0.5*I*p_chain->zeta())  , - p_chain->length());

					// product over other strings
					for (int k=0; k < p_base->numberTypes(); ++k)
					for (int beta=0; beta < p_base->numberStringsOfType(k); ++beta) {
						if ((k==j) && (alpha==beta)) continue;
						for (int b=0; b< p_chain->stringLength(k); ++b) {
							complex<REAL> root_k = root (k, beta, b);
							deviation[a] *= sinh (root_j - root_k + I* p_chain->zeta()) / sinh (root_j - root_k - I*p_chain->zeta());
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
	REAL XXZ_State::longitudinalFormFactor (const State& state_mu) const
	{
		// different subspaces do not overlap
		if (state_mu.p_base->numberDown() != p_base->numberDown()) return 0.0;
		return 0.25*p_chain->length() * exp(XXZ_State::lndetHm2P(state_mu));
	}



	/** transverse form factor |< GS | S- | state >|^2 **/
	REAL XXZ_State::transverseFormFactor (const State& state_mu) const
	{
		// different subspaces do not overlap
		if (state_mu.p_base->numberDown() != p_base->numberDown()+1) return 0.0;
		return p_chain->length() * exp(XXZ_State::lndetHminus(state_mu));
	}


	/** H minus 2 P determinant expression **/
	REAL XXZ_State::lndetHm2P (const State& state_mu) const
	{
		const char* here = "XXZ_State::lndetHm2P";
		if (state_mu.p_base->numberRoots() != p_base->numberRoots()) throw Exception (here, exc_NumberDown);

		// lambda, the right state of the matrix element
		int number_roots = p_base->numberRoots();
		complex<REAL>* root_lambda = new complex<REAL> [number_roots];
		complex<REAL>* sinh_lambda = new complex<REAL> [number_roots];
		complex<REAL>* cosh_lambda = new complex<REAL> [number_roots];
		for (int index=0, j=0; j< p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha)
		for (int a=0; a< p_chain->stringLength(j); ++a, ++index) {
			root_lambda[index] = root(j, alpha, a);
			sinh_lambda[index] = sinh(root_lambda[index]);
			cosh_lambda[index] = cosh(root_lambda[index]);
		}

		// mu, the left state of the matrix element
		complex<REAL>* root_mu = new complex<REAL> [number_roots];
		complex<REAL>* sinh_mu = new complex<REAL> [number_roots];
		complex<REAL>* cosh_mu = new complex<REAL> [number_roots];
		for (int index=0, j=0; j< state_mu.p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < state_mu.p_base->numberStringsOfType(j); ++alpha)
		for (int a=0; a< p_chain->stringLength(j); ++a, ++index) {
			root_mu[index] = state_mu.root(j, alpha, a);
			sinh_mu[index] = sinh(root_mu[index]);
			cosh_mu[index] = cosh(root_mu[index]);
		}
		// we only need these norm values at the very end, but calculate them now
		// so that there's no risk of having both the gaudin matrix
		// and hminus in memory at the same time.
		REAL norm_mu = state_mu.norm();
		REAL norm_lambda = norm();

		// these are all over the place.
		REAL cos_zeta = cos(p_chain->zeta());
		REAL sin_zeta = sin(p_chain->zeta());

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
				for (int k = 0; k < number_roots; ++k) {

					complex<REAL> sh_lk_m_lb = sinh_lambda[k] * cosh_lambda[b] - cosh_lambda[k] * sinh_lambda[b];
					complex<REAL> ch_lk_m_lb = cosh_lambda[k] * cosh_lambda[b] - sinh_lambda[k] * sinh_lambda[b];
					product_f0 *= sh_lk_m_lb * cos_zeta - I* ch_lk_m_lb * sin_zeta;
					if (k != b) product_f1 *= sh_lk_m_lb;
					if ((k != b + 1) || (b-csi == cs_length-1))
						product_f2 *= sh_lk_m_lb * cos_zeta + I*ch_lk_m_lb * sin_zeta;

					complex<REAL> sh_mk_m_lb = sinh_mu[k] * cosh_lambda[b] - cosh_mu[k] * sinh_lambda[b];
					complex<REAL> ch_mk_m_lb = cosh_mu[k] * cosh_lambda[b] - sinh_mu[k] * sinh_lambda[b];
					product_g0 *= sh_mk_m_lb * cos_zeta - I* ch_mk_m_lb * sin_zeta;
					product_g1 *= sh_mk_m_lb;
					product_g2 *= sh_mk_m_lb * cos_zeta + I* ch_mk_m_lb * sin_zeta;

					complex<REAL> sh_ma_m_mk = sinh_mu[b] * cosh_mu[k] - cosh_mu[b] * sinh_mu[k];
					complex<REAL> ch_ma_m_mk = cosh_mu[b] * cosh_mu[k] - sinh_mu[b] * sinh_mu[k];
					/**
					// filter out the zeroes in the denominator (need better implementation)
					// NOTE: watch out here! zero occurs for strings in the left (mu) state.
					// The only problem that arises is the zero in the denominator here;
					// however, in that case, it cancels against a 1/0 in the unreduced Gaudin matrix.
					// if we use the reduced Gaudin Matrix for the norm (as happes automatically when we feed a string state to the left)
					// we should leave this factor 0 out and get the correct finite result.
					**/
					if (::norm(sh_ma_m_mk*cos_zeta + I* ch_ma_m_mk*sin_zeta) > EQUAL_THRESHOLD)
						prefactor *= sh_ma_m_mk*cos_zeta + I* ch_ma_m_mk*sin_zeta;
					if (!(b+1==k) || (cs_length-1==vertical_b) )  {
						prefactor *= sh_lk_m_lb*cos_zeta + I*ch_lk_m_lb *sin_zeta;
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
					if (abs(root_mu[a] - root_lambda[b]) < EQUAL_THRESHOLD) throw Exception (here, exc_Equal);
					if (1==cs_length) {
						// 1+ and 1-
						complex<REAL> sh_ma_m_lcs = sinh_mu[a] * cosh_lambda[csi] - cosh_mu[a] * sinh_lambda[csi];
						complex<REAL> ch_ma_m_lcs = cosh_mu[a] * cosh_lambda[csi] - sinh_mu[a] * sinh_lambda[csi];
						matrix_hm2p[a][b] =
							1.0/ (  (sh_ma_m_lcs * cos_zeta - I* ch_ma_m_lcs * sin_zeta) * sh_ma_m_lcs  )
							+ exp(ln_g2[csi] - ln_f2[csi] - ln_g0[csi] + ln_f0[csi])
							  / ( sh_ma_m_lcs * (sh_ma_m_lcs * cos_zeta + I* ch_ma_m_lcs * sin_zeta))
							- 2.0 * exp(ln_f0[csi] - ln_g0[csi]) / (sq(sinh_mu[a]) + sq(sf_half_n_anis[1]));
					}
					else if (vertical_b < cs_length-1) {
						// strings, element for "b < n"
						matrix_hm2p[a][b] = exp(-ln_g0[b])/( sinh(root_mu[a] - root_lambda[b])
													* sinh(root_mu[a] - root_lambda[b] + I*p_chain->zeta()) );
					}
					else {
						// strings, element for "b==n"
						// in these comments, subscripts indicate index as in paper; brackets indicate subscript as in this code.
						// lambda_1 == current_string_lambda[0]; I don't use lambda_0 and lambda_{n+1} since they're not real rapidities
						complex<REAL> sum_h (0.0, 0.0), sum_p (0.0, 0.0);

						// j=0 term ( a non-rapidity; remember, lambda_j with j==1 refers to the current root_lambda[b] )
						// paper: K * ( G_0 G_1 ) / ( F_0 F_1 )
						complex<REAL> sh_ma_m_lcs = sinh_mu[a] * cosh_lambda[csi] - cosh_mu[a] * sinh_lambda[csi];
						complex<REAL> ch_ma_m_lcs = cosh_mu[a] * cosh_lambda[csi] - sinh_mu[a] * sinh_lambda[csi];
						sum_h += exp(ln_g0[csi] + ln_g1[csi] - ln_f0[csi] - ln_f1[csi]
									- log( ( sh_ma_m_lcs * cos_zeta - I* ch_ma_m_lcs * sin_zeta ) * sh_ma_m_lcs )
									);

						// 0 < j < n terms: -d/d lambda K_j * ( G_j G_{j+1} ) / ( F_j F_{j+1} )
						for (int s_i = 0; s_i < cs_length-1; ++s_i) {
							complex<REAL> sh_ma_m_lbpi = sinh_mu[a] * cosh_lambda[csi+s_i] - cosh_mu[a] * sinh_lambda[csi+s_i];
							complex<REAL> ch_ma_m_lbpi = cosh_mu[a] * cosh_lambda[csi+s_i] - sinh_mu[a] * sinh_lambda[csi+s_i];
							sum_p += exp( ln_g1[csi+s_i] - ln_f1[csi+s_i] );
							sum_h -= exp( ln_g1[csi+s_i] + ln_g2[csi+s_i] - ln_f1[csi+s_i] - ln_f2[csi+s_i]
										+ log(  ( ch_ma_m_lbpi / sh_ma_m_lbpi
												  +	(ch_ma_m_lbpi * cos_zeta + I* sh_ma_m_lbpi * sin_zeta)
													/ (sh_ma_m_lbpi * cos_zeta + I* ch_ma_m_lbpi * sin_zeta) )
												/ ( sh_ma_m_lbpi * (sh_ma_m_lbpi * cos_zeta + I* ch_ma_m_lbpi * sin_zeta)  ) )  );
						}

						// j=n term : K * ( G_n G_{n+1} ) / ( F_n F_{n+1} )
						complex<REAL> sh_ma_m_lbpi = sinh_mu[a] * cosh_lambda[csi+cs_length-1] - cosh_mu[a] * sinh_lambda[csi+cs_length-1];
						complex<REAL> ch_ma_m_lbpi = cosh_mu[a] * cosh_lambda[csi+cs_length-1] - sinh_mu[a] * sinh_lambda[csi+cs_length-1];
						sum_p += exp( ln_g1[csi+cs_length-1] - ln_f1[csi+cs_length-1] );
						sum_h += exp( ln_g1[csi+cs_length-1] + ln_g2[csi+cs_length-1] - ln_f1[csi+cs_length-1] - ln_f2[csi+cs_length-1]
										- log( sh_ma_m_lbpi * (sh_ma_m_lbpi * cos_zeta + I* ch_ma_m_lbpi * sin_zeta)) 	);

						complex<REAL> ln_normalisation = ln_f0[csi] + ln_f1[csi] - ln_g0[b];	// N *= F1/G1 (now it is F0*F1/G1)
						for (int s_i = 1; s_i < cs_length; ++s_i)
							ln_normalisation += ln_g1[csi+s_i];
						ln_normalisation -= ln_g1[csi+cs_length-1];

						sum_p /= sq(sinh_mu[a]) + sq(sf_half_n_anis[1]);
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

		REAL lndet_hm2p = lndetDestroy (matrix_hm2p);
		// apply the a-independent factors (i.e. we have divided the matrix columns by factors which we now apply to the determinant)
		for (int b=0; b < number_roots; ++b)
			lndet_hm2p += real(ln_g0[b]);  // ln(abs(z)) == real(ln(z))

		REAL ln_form_factor = 0.0;
		for (int index=0; index < number_roots; ++index)
			ln_form_factor += 2.0 * log(abs(sinh(root_mu[index] - 0.5*I*p_chain->zeta()))) - 2.0 * log(abs(sinh(root_lambda[index] - 0.5*I*p_chain->zeta())));
		ln_form_factor += 2.0 * number_roots * log( sin(p_chain->zeta()) );
		ln_form_factor += 2.0 * lndet_hm2p + ln_pre - norm_mu - norm_lambda;

		delete[] root_lambda; delete[] sinh_lambda; delete[] cosh_lambda;
		delete[] root_mu; delete[] sinh_mu; delete[] cosh_mu;
		delete[] ln_f0; delete[] ln_f1; delete[] ln_f2;
		delete[] ln_g0; delete[] ln_g1; delete[] ln_g2;
		return ln_form_factor;
	}


	/** H-minus determinant expression **/
	REAL XXZ_State::lndetHminus (const State& state_mu) const
	{
	 	const char* here = "XXZ_State::lndetHminus";
		if (state_mu.p_base->numberRoots() != p_base->numberRoots()+1) throw Exception (here, exc_NumberDown);

		// lambda, the right state of the matrix element
		int number_roots = p_base->numberRoots();
		complex<REAL>* root_lambda = new complex<REAL> [number_roots];
		complex<REAL>* sinh_lambda = new complex<REAL> [number_roots];
		complex<REAL>* cosh_lambda = new complex<REAL> [number_roots];
		for (int index=0, j=0; j< p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < p_base->numberStringsOfType(j); ++alpha)
		for (int a=0; a< p_chain->stringLength(j); ++a, ++index) {
			root_lambda[index] = root(j, alpha,a);
			sinh_lambda[index] = sinh(root_lambda[index]);
			cosh_lambda[index] = cosh(root_lambda[index]);
		}

		// mu, the left state of the matrix element
		int number_roots_mu = state_mu.p_base->numberRoots();
		complex<REAL>* root_mu = new complex<REAL> [number_roots_mu];
		complex<REAL>* sinh_mu = new complex<REAL> [number_roots_mu];
		complex<REAL>* cosh_mu = new complex<REAL> [number_roots_mu];
		for (int index=0, j=0; j< state_mu.p_base->numberTypes(); ++j)
		for (int alpha=0; alpha < state_mu.p_base->numberStringsOfType(j); ++alpha)
		for (int a=0; a< p_chain->stringLength(j); ++a, ++index) {
			root_mu[index] = state_mu.root(j, alpha,a);
			sinh_mu[index] = sinh(root_mu[index]);
			cosh_mu[index] = cosh(root_mu[index]);
		}

		// we only need these norm values at the very end, but calculate them now
		// so that there's no risk of having both the gaudin matrix
		// and hminus in memory at the same time.
		REAL norm_mu = state_mu.norm();
		REAL norm_lambda = norm();

		// these are all over the place.
		REAL cos_zeta = cos(p_chain->zeta());
		REAL sin_zeta = sin(p_chain->zeta());

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
						complex<REAL> sh_lk_m_lb = sinh_lambda[k] * cosh_lambda[b] - cosh_lambda[k] * sinh_lambda[b];
						complex<REAL> ch_lk_m_lb = cosh_lambda[k] * cosh_lambda[b] - sinh_lambda[k] * sinh_lambda[b];
						product_f0 *= sh_lk_m_lb * cos_zeta - I* ch_lk_m_lb * sin_zeta;
						if (k != b) product_f1 *= sh_lk_m_lb;
						if ((k != b + 1) || (b-csi == cs_length-1))
							product_f2 *= sh_lk_m_lb * cos_zeta + I*ch_lk_m_lb * sin_zeta;
						if (!(b+1==k) || (cs_length-1==vertical_b) )  {
							prefactor *= sh_lk_m_lb*cos_zeta + I*ch_lk_m_lb *sin_zeta;
						}
					}
					complex<REAL> sh_mk_m_lb = sinh_mu[k] * cosh_lambda[b] - cosh_mu[k] * sinh_lambda[b];
					complex<REAL> ch_mk_m_lb = cosh_mu[k] * cosh_lambda[b] - sinh_mu[k] * sinh_lambda[b];
					product_g0 *= sh_mk_m_lb * cos_zeta - I* ch_mk_m_lb * sin_zeta;
					product_g1 *= sh_mk_m_lb;
					product_g2 *= sh_mk_m_lb * cos_zeta + I* ch_mk_m_lb * sin_zeta;

					complex<REAL> sh_ma_m_mk = sinh_mu[b] * cosh_mu[k] - cosh_mu[b] * sinh_mu[k];
					complex<REAL> ch_ma_m_mk = cosh_mu[b] * cosh_mu[k] - sinh_mu[b] * sinh_mu[k];
					// filter out the zeroes in teh denominator (need better implementation)
					if (::norm(sh_ma_m_mk*cos_zeta + I* ch_ma_m_mk*sin_zeta) > EQUAL_THRESHOLD)
						prefactor *= sh_ma_m_mk*cos_zeta + I* ch_ma_m_mk*sin_zeta;

					// gather results in log (which avoids overflows)
					// but not too often as logs take time.
					if (!(k%50) || (k==number_roots)) {
						ln_f0[b]+= log(product_f0);
						ln_f1[b]+= log(product_f1);
						ln_f2[b]+= log(product_f2);
						ln_g0[b]+= log(product_g0);
						ln_g1[b]+= log(product_g1);
						ln_g2[b]+= log(product_g2);
						ln_pre -= real(log( prefactor ));
						product_f0 = product_f1 = product_f2 = 1.0;
						product_g0 = product_g1 = product_g2 = 1.0;
						prefactor=1.0;
					}
				}
			}
		}

		// prefactor for last column b == M+1
		complex<REAL> prefactor = 1.0;
		for (int k = 0; k <= number_roots; ++k) {
			complex<REAL> sh_ma_m_mk = sinh_mu[number_roots_mu-1] * cosh_mu[k] - cosh_mu[number_roots_mu-1] * sinh_mu[k];
			complex<REAL> ch_ma_m_mk = cosh_mu[number_roots_mu-1] * cosh_mu[k] - sinh_mu[number_roots_mu-1] * sinh_mu[k];
			if (::norm(sh_ma_m_mk*cos_zeta + I* ch_ma_m_mk*sin_zeta) > EQUAL_THRESHOLD)
				prefactor *= sh_ma_m_mk*cos_zeta + I* ch_ma_m_mk*sin_zeta;
			// gather results in log (which avoids overflows)
			// but not too often as logs take time.
			if (!(k%50) || (k==number_roots)) {
				ln_pre -= real(log( prefactor ));
				prefactor=1.0;
			}
		}

		Square< complex<REAL> > matrix_hminus (number_roots_mu);
		// type_a runs over all string types
		// alpha runs over all instances of a given string length
		// vertical_a runs over all roots in a given string instance.
		// a runs over all roots (first/row matrix index, goes with mu)
		for (int a=0, type_a=0; type_a< state_mu.p_base->numberTypes(); ++type_a)
		for (int alpha=0; alpha < state_mu.p_base->numberStringsOfType(type_a); ++alpha)
		for (int vertical_a=0; vertical_a < state_mu.p_chain->stringLength(type_a); ++vertical_a, ++a) {

			// last column of H-
			matrix_hminus[a][number_roots_mu-1] =	1.0 / (sq(sinh_mu[a]) + sq(sf_half_n_anis[1]));

			// indices type_b, beta, vertical_B, b (second/column matrix index, goes with lambda): analogous to a.
			for (int b=0, type_b=0; type_b<p_base->numberTypes(); ++type_b)
			for (int beta=0; beta < p_base->numberStringsOfType(type_b); ++beta) {
				int csi = b; 									// current string index: always points to the first root in the current string
				int cs_length = p_chain->stringLength(type_b); 	// length of the current string
				for (int vertical_b = 0; vertical_b < cs_length; ++vertical_b, ++b) {
					if (abs(root_mu[a] - root_lambda[b]) < EQUAL_THRESHOLD) throw Exception (here, exc_Equal);
					if (1==cs_length) {
						// 1+ and 1- rapidities
						complex<REAL> sh_ma_m_lcs = sinh_mu[a] * cosh_lambda[csi] - cosh_mu[a] * sinh_lambda[csi];
						complex<REAL> ch_ma_m_lcs = cosh_mu[a] * cosh_lambda[csi] - sinh_mu[a] * sinh_lambda[csi];
						matrix_hminus[a][b] =
							1.0/ ( ( sh_ma_m_lcs * cos_zeta - I* ch_ma_m_lcs * sin_zeta ) * sh_ma_m_lcs )
							+ exp(- ln_g0[csi]  + ln_g2[csi] - ln_f2[csi] + ln_f0[csi] )
								/ ( sh_ma_m_lcs * (sh_ma_m_lcs * cos_zeta + I* ch_ma_m_lcs * sin_zeta))  ;
					}
					else if (vertical_b < cs_length-1) {
						// strings, element for "b<n"
						matrix_hminus[a][b] = exp(-ln_g0[b])/( sinh(root_mu[a] - root_lambda[b]) * sinh(root_mu[a] - root_lambda[b] + I*p_chain->zeta()) );
					}
					else {
						// strings, element for "b==n"
						// in these comments, subscripts indicate index as in lateX notes; brackets indicate subscript as in this code.
						// lambda_1 == current_string_lambda[0] ; I don't use lambda_0 and lambda_{n+1} since they're not real rapidities
						complex<REAL> sum_h (0.0, 0.0);

						// j=0 term ( a non-rapidity; remember, lambda_j with j==1 refers to the current root_lambda[b] )
						// paper: K * ( G_0 G_1 ) / ( F_0 F_1 )
						complex<REAL> sh_ma_m_lcs = sinh_mu[a] * cosh_lambda[csi] - cosh_mu[a] * sinh_lambda[csi];
						complex<REAL> ch_ma_m_lcs = cosh_mu[a] * cosh_lambda[csi] - sinh_mu[a] * sinh_lambda[csi];
						sum_h += 	exp(ln_g0[csi] + ln_g1[csi] - ln_f0[csi] - ln_f1[csi]
										- log( ( sh_ma_m_lcs * cos_zeta - I* ch_ma_m_lcs * sin_zeta ) * sh_ma_m_lcs ) );

						// 0 < j < n terms: -d/d lambda K_j * ( G_j G_{j+1} ) / ( F_j F_{j+1} )
						for (int s_i = 0; s_i < cs_length-1; ++s_i) {
							// sinh and cosh of root_mu[a] - root_lambda[b]
							complex<REAL> sh_ma_m_lbpi = sinh_mu[a] * cosh_lambda[csi+s_i] - cosh_mu[a] * sinh_lambda[csi+s_i];
							complex<REAL> ch_ma_m_lbpi = cosh_mu[a] * cosh_lambda[csi+s_i] - sinh_mu[a] * sinh_lambda[csi+s_i];
							sum_h -= exp( ln_g1[csi+s_i] + ln_g2[csi+s_i] - ln_f1[csi+s_i] - ln_f2[csi+s_i]
										+ log( ( ch_ma_m_lbpi / sh_ma_m_lbpi
													+ 	(ch_ma_m_lbpi * cos_zeta + I* sh_ma_m_lbpi * sin_zeta)
													    / (sh_ma_m_lbpi * cos_zeta + I* ch_ma_m_lbpi * sin_zeta)	)
												/ ( sh_ma_m_lbpi * (sh_ma_m_lbpi * cos_zeta + I* ch_ma_m_lbpi * sin_zeta)  )		) );
						}

						// j=n term : K * ( G_n G_{n+1} ) / ( F_n F_{n+1} )
						complex<REAL> sh_ma_m_lbpi = sinh_mu[a] * cosh_lambda[csi+cs_length-1] - cosh_mu[a] * sinh_lambda[csi+cs_length-1];
						complex<REAL> ch_ma_m_lbpi = cosh_mu[a] * cosh_lambda[csi+cs_length-1] - sinh_mu[a] * sinh_lambda[csi+cs_length-1];
						sum_h += exp( ln_g1[csi+cs_length-1] + ln_g2[csi+cs_length-1] - ln_f1[csi+cs_length-1] - ln_f2[csi+cs_length-1]
												- log( sh_ma_m_lbpi * (sh_ma_m_lbpi * cos_zeta + I* ch_ma_m_lbpi * sin_zeta)) );

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

		REAL lndet_hminus = lndetDestroy (matrix_hminus);
		for (int b=0; b < number_roots; ++b)
			lndet_hminus += real(ln_g0[b]);  // ln(abs(z)) == real(ln(z))

		REAL ln_form_factor = 0.0;
		for (int index=0; index < number_roots; ++index)
			ln_form_factor += 2.0 * log(abs(sinh(root_mu[index] - 0.5*I*p_chain->zeta()))) - 2.0 * log(abs(sinh(root_lambda[index] - 0.5*I*p_chain->zeta())));
		ln_form_factor += 2.0 * log(abs(sinh(root_mu[number_roots_mu-1] - 0.5*I*p_chain->zeta())));
		ln_form_factor += (2.0 * number_roots + 2.0) * log( sin(p_chain->zeta()) );
		ln_form_factor += 2.0 * lndet_hminus + ln_pre - norm_mu - norm_lambda;

		delete[] root_lambda; delete[] sinh_lambda; delete[] cosh_lambda;
		delete[] root_mu; delete[] sinh_mu; delete[] cosh_mu;
		delete[] ln_f0; delete[] ln_f1; delete[] ln_f2;
		delete[] ln_g0; delete[] ln_g1; delete[] ln_g2;
		return ln_form_factor;
	}


