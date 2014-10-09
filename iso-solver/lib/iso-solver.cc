#include "iso-solver.h"


/** exceptions **/

const string exc_WrongPairWidth = "narrow pair mistaken for wide or vice versa";


	
	
/** default deviated solution policies **/

// maximum string deviation we accept before trying to solve for deviance
double IsoSolver::deviation_threshold = 1e-20;

// maximum number of iteratons in deviance solution
int IsoSolver::max_iter_deviance = 1000;
// precision to be reached in deviance iteration
double IsoSolver::deviance_precision = 1e-28; 
// how bad 'convergence' may become before we call it quits
double IsoSolver::solve_bailout = 1e+6;
// damping for deviance	
double IsoSolver::damping_deviance = 0.5;
// damping for solveSymmetric()
double IsoSolver::damping_symmetric = 0.5;



const double IsoSolver::no_convergence = sq(BOOMBANOOMBA);
const double IsoSolver::not_calculated = -sq(BOOMBANOOMBA);

const double IsoSolver::start_deviance = 1e-20;
// epsilon in the check for deviance=-0.5. NOTE that this need not be very small, as we don't really get 'near-collapses'	
const double IsoSolver::epsilon_collapse = 1e-5;
// below which the norm would give a non-finite result if we would use the complex roots (i.e. we should use the string equations for the norm)
const double IsoSolver::epsilon_deviation_norm = 1e-16;
// how close roots may come to i/2 before we consider it to be exactly that for the calculation of energies
const double IsoSolver::epsilon_energy = 1e-10;

	



/** xi: argument-like function or 'continuous arctan', see notes. **/

// argument-like function or 'continuous arctan', see notes.
inline long double xi (const long double epsilon, const long double delta) 
{
	return (atan(epsilon/delta) + ((delta<0.0)?PI*sgn(epsilon):0));
}	

	
/** IsoSolver **/	



    IsoSolver::IsoSolver (const int chain_length, const IsoConfiguration& on_base, const Strip<int>& quantum)
        :   IsoBetheTakahashiSolver(chain_length, on_base, quantum), 
            deviance(on_base, start_deviance), aberration(on_base, start_deviance),
            hold(on_base, false), hold_deviation(on_base, false)
    { 	 
        resetDeviances();
    }
	
	
	/* set free rapidies */
	void IsoSolver::reset()
	{
		IsoBetheTakahashiSolver::reset();
        resetDeviances();
	}


	void IsoSolver::resetDeviances()
	{
		    for (int j=0; j < config_.numberTypes(); ++j) {  
		    // set sensible initial conditions
		    int length_j = config_.stringLength(j);
		    for (int alpha=0; alpha < config_.numberStringsOfType(j); ++alpha) {
			    if (length_j%2) {
				    deviance (j, alpha, (length_j-1)/2) = 0;
				    aberration (j, alpha, (length_j-1)/2) = 0;
			    }
			    else {
				    aberration (j, alpha, length_j/2-1) = 0;
				    aberration (j, alpha, length_j/2) = 0;
			    }
			    for (int a=(config_.stringLength(j)+1)/2; a < config_.stringLength(j); ++a) {
				    deviance(j, alpha, a) = -start_deviance;
			    }
		    }
	    }
	
	}
	
	
    /* roots of the bethe equation */
    vector< complex<double> > IsoSolver::roots (void) const
    {
	    vector< complex<double> > root (config_.numberRoots());
	    for (int index=0, j=0; j < config_.numberTypes(); ++j)  
	    for (int alpha=0; alpha < config_.numberStringsOfType(j); ++alpha) 
	    for (int a=0; a < config_.stringLength(j); ++a, ++index) {
		    root[index] =  IsoSolver::root(j, alpha, a);
	    }
	    return root;
    }


    /* roots of the bethe equation */
    complex<double> IsoSolver::root (const int j, const int alpha, const int a) const
    {	
	    return rapidity(j, alpha) + aberration(j,alpha,a) 	+  ( 0.5 *(config_.stringLength(j)-1-2*a) + deviance(j,alpha,a) )*I;	
    }



	/** solve by iteration **/
	bool IsoSolver::solve(void) 
	{
		const char* here = "IsoSolver::solve";
		
		// first solve the equations without deviance
		if (!IsoBetheTakahashiSolver::solve()) return false;	

cerr<<"bt roots "<< roots() <<endl;	
for (int j=0; j<config_.numberTypes(); ++j)
for (int alpha=0; alpha<config_.numberStringsOfType(j); ++alpha) {
    cerr<<sumOfBethe2xJ(j,alpha)<<";";
}
                
		// low deviation? don't look any further
		if (IsoBetheTakahashiSolver::stringDeviation() < deviation_threshold) return true;
	
		// reset convergence to indicate we're none of your converged states
		convergence = no_convergence;

		// odd symmetric or inadmissible states: use different logic
		if (!admissible())  return solveSingular();
		else if (oddSymmetric()) return solveSymmetric();
		
/*** RLH2013 ***/
        else {
        
            //find apriori colliding roots 
            for (int j=0; j < config_.numberTypes(); ++j) 
		    for (int alpha=0; alpha < config_.numberStringsOfType(j); ++alpha) {
		        // only check central roots in odd strings
		        if (0==config_.stringLength(j)%2) continue;
		        
		        for (int k=j+1; k < config_.numberTypes(); ++k) 
		        for (int beta=0; beta < config_.numberStringsOfType(k); ++beta) {
		            if (0==config_.stringLength(k)%2) continue;
		        
		            if (abs( rapidity(j,alpha) - rapidity(k,beta)) < 1e-5)  {
cerr<<" apriori bethe-takahashi collision "<<endl;
                        // bump
		                rapidity(0,0 ) += 20;
		cerr<<roots()<<endl;
		                // can we use symmetric-state solver for this kind?
		            }
		        }
		    }            
                
        }
/*** /RLH2013 ***/
        
        
		// iterate to convergence
		while (abs(convergence) > precision && iterations < max_iterations) {
			// our iterate method is aware of deviance and increments iterations
 			iterate();
 			// bail out if we're diverging too much
 			if (abs(convergence) > solve_bailout) return false;
		}
		if (abs(convergence) > precision) return false;

		// calculate deviation te see if we have an actual solution or a LHS==-RHS type result
		//if (IsoSolver::stringDeviation() > DEVIANCE_PRECISION) return false;
		// if the deviance equation and the Bethe equation both converge, we should have a valid solution if delta is reasonable

		for (int j=0; j < config_.numberTypes(); ++j) 
		for (int alpha=0; alpha < config_.numberStringsOfType(j); ++alpha) 
		for (int a=0; a < config_.stringLength(j)-1; ++a) 
		{
			// check for collapses on a string.
			// note that odd symmetric complexes have very close real central roots, but are okay.
			// since they're on different strings, they won't trigger this test.
			if (abs(root(j, alpha, a) - root(j, alpha, a+1)) < epsilon_collapse) 
			{
				// collapse signals extra real solutions. new approach.
				
/** RLH2008 **/	
cerr<<"collapse "<<j<<' '<<alpha<<' '<<a<<' '<<deviance(j, alpha,a)<<' '<<aberration(j,alpha,a)<<endl;					
cerr<<roots();
				// we assume a two-string. longer collapsing strings not yet supported.
				hold(j, alpha) = true;
				hold_deviation(j,alpha,a) = true;
				hold_deviation(j,alpha,a+1) = true;
				
				convergence = no_convergence;
				while (abs(convergence) > precision && iterations < max_iterations) {
					// our iterate method is aware of deviance and increments iterations

					if (!findExtraReal(j, alpha, a)) return false;
					iterate();
					if (abs(convergence) > solve_bailout) return false;
				}
				if (abs(convergence) > precision) return false;
/** /RLH2008 **/				
				
				
			}
		}
		return true;
	}
		

	/** calculate a new value for the deviance (delta) and aberration (epsilon) **/
	bool IsoSolver::getNewDeviances (FullStrip<double>& new_deviance, FullStrip<double>& new_aberration, const int j, const int alpha) const
	{
		
		const string here = "IsoSolver::getNewDeviances";
		
		// NOTE that int a in program runs from 0..len_j-1 , in notes from 1..n_j

		// the real rapidities: they deviate not, nor do they aberr.
		int length_j = config_.stringLength(j);
 		if (1==length_j) 		{
			new_deviance(j,alpha,0) = new_aberration(j,alpha,0) = 0;
			return true;
		}
		double lam = rapidity(j, alpha);
		double big_n = chain_length_;

		long double theta = 0.0;
		long double log_r = 0.0;
		long double del_term[length_j/2];
		long double eps_term[length_j/2];
		for (int a=0; a<length_j/2; ++a) {
			double del = deviance (j, alpha, a); 
			double eps = aberration (j, alpha, a);
		
			// contribution M mod 2 from sum of conjugate quantum numbers
			theta -= PI*config_.numberRoots();		
			
			// kinetic theta
			theta += big_n * ( 
					  xi(2.0*(lam+eps), length_j - 2.0*a + 2.0*del) 
					+ xi(2.0*(lam+eps), -length_j + 2.0*(a+1) - 2.0*del)
				);
			log_r += big_n * (
 					   log(sq(lam+eps) + sq(0.5*(length_j) - (a+1) + del))
 					 - log(sq(lam+eps) + sq(0.5*(length_j) - (a) + del))
 				);

			for (int k=0; k < config_.numberTypes(); ++k) 
			for (int beta=0; beta < config_.numberStringsOfType(k); ++beta) {
				if (k==j && beta==alpha) {
					// theta_self minus dangerous terms		
					for (int b=0; b < length_j; ++b) {
						if (b==a) continue;
						double del_b = deviance (j, alpha, b); 
						double eps_b = aberration (j, alpha, b);
						if (b==a-1) {
							theta -= xi( (eps-eps_b), (2.0-del+del_b) );	
							log_r -= log( sq(eps-eps_b) + sq(2.0-del+del_b) );
						}
						else if (b==a+1) {
 							theta -= xi( (eps-eps_b), (2.0+del-del_b) );	
							log_r += log( sq(eps-eps_b) + sq(2.0+del-del_b) );
						}
						else {
							theta -=  xi( (eps-eps_b), (1.0-(a-b)+(del-del_b)) )
									+ xi( (eps-eps_b), (1.0+(a-b)-(del-del_b)) );

							log_r +=  log( 	sq(eps-eps_b)  +  sq(1.0-(a-b)+(del-del_b)) 	)
									- log( 	sq(eps-eps_b)  +  sq(1.0+(a-b)-(del-del_b)) 	);
						}
					}		
				}
				else {
					// theta_other
					int length_k = config_.stringLength(k);
					double lam_k = rapidity(k, beta);
					for (int b=0; b < length_k; ++b) {
						double del_k = deviance (k, beta, b); 
						double eps_k = aberration (k, beta, b);
					
						theta -=  xi( lam-lam_k + eps-eps_k, 1.0 + 0.5*(length_j-length_k) - (a-b) + (del-del_k) )
								+ xi( lam-lam_k + eps-eps_k, 1.0 - 0.5*(length_j-length_k) + (a-b) - (del-del_k) );
							
						log_r +=  log( sq(lam-lam_k+eps-eps_k)  +  sq( 1.0+0.5*(length_j-length_k) - (a-b) + (del-del_k)) )
								- log( sq(lam-lam_k+eps-eps_k)  +  sq(-1.0+0.5*(length_j-length_k) - (a-b) + (del-del_k)) );
					}
				}
			}
			// theta is now sum b=1..a theta[b]
			// log_r is now sum b=1..a log_r[b]
			del_term[a] = -cos(theta)*exp(0.5*log_r);
			eps_term[a] = sin(theta)*exp(0.5*log_r);
		}		
		if (length_j%2) {
			// odd n_j
			for (int a=0; a<(length_j-1)/2; ++a) {
				for (int b=a; b<length_j/2; ++b) {
					new_deviance (j,alpha,a) += del_term[b];
					new_aberration (j,alpha,a) += eps_term[b];
				}
			}
			// string center has no deviations.
			new_deviance(j,alpha, (length_j-1)/2) = 0.0;
			new_aberration(j, alpha, (length_j-1)/2) = 0.0;
			// set the conjugates
			for (int a=(length_j+1)/2; a<length_j;++a) {
				new_deviance(j,alpha,a) = - new_deviance(j,alpha,length_j-a-1);
				new_aberration(j,alpha,a) = new_aberration(j,alpha,length_j-a-1);
			}
		}
		else {
			// even n_j. 
			// this is more involved, see notes
			
			// first, the regular terms.
			for (int a=0; a<length_j/2-1; ++a) {
				for (int b=a; b<length_j/2-1; ++b) {
					new_deviance (j,alpha,a) += del_term[b];
					new_aberration (j,alpha,a) += eps_term[b];
				}
			}
			// determine width of inner pair
			int sum_sign = 0;
			for (int k=0; k < config_.numberTypes(); ++k) 
			for (int beta=0; beta < config_.numberStringsOfType(k); ++beta) {
				if (j==k && alpha==beta) continue;
				int length_k = config_.stringLength(k);
				
				sum_sign += isgn(rapidity(j, alpha)-rapidity(k, beta))*(length_j*length_k - 2*min(length_j, length_k) + (j==k));
			}
 			int sign_inner = ((quantum_number(j, alpha) + chain_length_*(length_j - 1) - length_j*(config_.numberRoots()+2) - 2 - sum_sign)/2)%2 	? 1 : -1;
			// apply the inner-sign term.
			for (int a=0; a<length_j/2; ++a) {
				// log_r should be sum b=1...n/2 log_r[b] now
				new_deviance(j,alpha, a) += 0.5*sign_inner*exp(0.5*log_r); 
			}

			
			// inner pair has no aberration
			new_aberration(j, alpha, length_j/2-1)=0;
			new_aberration(j, alpha, length_j/2)=0;
			new_deviance(j, alpha, length_j/2)=-new_deviance(j, alpha, length_j/2-1);
			
			// set the conjugates
			for (int a=length_j/2+1; a<length_j;++a) {
				new_deviance(j,alpha,a) = - new_deviance(j,alpha,length_j-a-1);
				new_aberration(j,alpha,a) = new_aberration(j,alpha,length_j-a-1);
			}
		}
		
		for (int a=0; a<length_j; ++a) 
			new_deviance(j,alpha, a) = (1.0-damping_deviance)*new_deviance(j,alpha, a) + damping_deviance*deviance(j,alpha,a);
		
		// enforce holds
		for (int a=0; a<length_j;++a) {
			if (hold_deviation(j, alpha, a)) {
				new_deviance(j,alpha,a) = deviance(j, alpha, a);
				new_aberration(j,alpha,a) = aberration(j,alpha, a);
			}
		}
		return true;
	}

	/** calculate a new value for the string centre (lambda) **/
	long double IsoSolver::newRapidity  (const int j, const int alpha) const
	{
		const char* here ="IsoSolver::newRapidity";
		int length_j = config_.stringLength(j);

		double lam = rapidity(j, alpha);
		long double scattering_term = 0.0;		
		for (int a=0; a<length_j; ++a) {
			double del = deviance (j, alpha, a);
			double eps = aberration (j, alpha, a);
		
			for (int k=0; k < config_.numberTypes(); ++k) {
				int length_k = config_.stringLength(k);
				for (int beta=0; beta < config_.numberStringsOfType(k); ++beta) {
					if ((j==k)&&(alpha==beta)) continue; 
					long double lj_m_lk_p_e = lam - rapidity(k, beta) + eps;	
					if (1==length_k) 
						scattering_term += 	  xi( lj_m_lk_p_e,  1.0 - (0.5*(length_j-1) - a + del) )
											+ xi( lj_m_lk_p_e,  1.0 + (0.5*(length_j-1) - a + del) );
					else {
						scattering_term +=    xi( lj_m_lk_p_e - aberration(k,beta,0),  1.0 - (0.5*(length_j-length_k) + 0 - a + del - deviance(k,beta,0)  ) )
					 						+ xi( lj_m_lk_p_e - aberration(k,beta,1),  1.0 - (0.5*(length_j-length_k) + 1 - a + del - deviance(k,beta,1)  ) )
					 						+ xi( 	lj_m_lk_p_e - aberration(k,beta,length_k-2),
					 								1.0 + (0.5*(length_j-length_k) + (length_k-2) - a + del - deviance(k,beta,length_k-2)  ) )
					 						+ xi( 	lj_m_lk_p_e - aberration(k,beta,length_k-1),
					 								1.0 + (0.5*(length_j-length_k) + (length_k-1) - a + del - deviance(k,beta,length_k-1)  ) );
					
						for (int b=0; b< length_k-2; ++b) 
							scattering_term +=    xi(  lj_m_lk_p_e - aberration(k,beta,b  ),  1.0 + ( 0.5*(length_j-length_k) + b   - a + del - deviance(k,beta,b  ) )  )	
												+ xi(  lj_m_lk_p_e - aberration(k,beta,b+2),  1.0 - ( 0.5*(length_j-length_k) + b+2 - a + del - deviance(k,beta,b+2) )  );
					}
				}
			}
		}
		// extract sum of Bethe J's from Takahashi I
		int sum_2x_bethe_quantum = quantum_number(j, alpha) + chain_length_*(length_j-1);
		for (int k=0; k < config_.numberTypes(); ++k) {
			int length_k = config_.stringLength(k);
			for (int beta=0; beta < config_.numberStringsOfType(k); ++beta) {
				if (j==k && beta==alpha) continue;
				sum_2x_bethe_quantum -= isgn(lam - rapidity(k, beta))*(length_j*length_k - 2*min(length_j, length_k) + (j==k));
			}
		}
		
		// 2\pi \sum_{a=1}^n J_a + \theta\sub{other}^a = \sum_{a=1}^n \theta\sub{kin}^a 
		long double phase = 0.5*(PI*sum_2x_bethe_quantum + scattering_term)/double(chain_length_);
		
		
		if (length_j%2) {
			// odd length
			long double sum_kin = 0.0;
			for (int a=0; a< (length_j-1)/2; ++a) {
				sum_kin += 	  xi(  (lam + aberration(j,alpha,a+1)),  ( 0.5*length_j - (a+1) + deviance(j,alpha,a+1))  )
							+ xi(  (lam + aberration(j,alpha,a)),  (-0.5*length_j + (a+1) - deviance(j,alpha,a))  );
			}
			return - aberration(j,alpha,0) + ( 0.5*length_j + deviance(j,alpha,0) )*tan( phase-sum_kin ) ;
		}
		else {
			// even length
			long double sum_kin = 0.0;
			for (int a=0; a< length_j/2-1; ++a) {
				sum_kin += 	  xi(  (lam + aberration(j,alpha,a+1)),  ( 0.5*length_j-(a+1) + deviance(j,alpha,a+1)) )
							+ xi(  (lam + aberration(j,alpha,a)),  (-0.5*length_j+(a+1) - deviance(j,alpha,a)) );
			}
			double tan_theta =  tan(  phase - sum_kin  ) ;
			double eps_1 = aberration(j, alpha, 0);
			double del_1 = deviance(j, alpha, 0);
			double del_n2 = deviance(j, alpha, length_j/2-1);
			double b = (eps_1*tan_theta+del_1-del_n2+0.5*length_j);
			double c = del_n2*(tan_theta*(0.5*length_j+del_1)-eps_1);
			// keep the sign as it was
			int sign_sqrt = isgn(2.0*lam*tan_theta+b);
			return (-b + sign_sqrt*sqrt( sq(b) - 4.0*tan_theta*c ))/(2.0*tan_theta);
		}

		return not_calculated;
	}
	
	
	
	
	/** one iteration of Bethe's equation **/
	void IsoSolver::iterate (void) 
	{
		string here = "IsoSolver::iterate"; 
		Strip<double> new_rapidity (config_);
		FullStrip<double> new_deviance (config_);
		FullStrip<double> new_aberration (config_);
		double square_diff = 0.0;	
		double one_over_n = 1.0/(1.0*chain_length_);
		long double new_rap = 0.0;
		
		for (int j = 0; j < config_.numberTypes(); ++j) 
		for (int alpha = 0; alpha < config_.numberStringsOfType(j); ++alpha) {
			if (!hold(j, alpha)) {
				// this should be finite
				double current_rapidity = rapidity(j, alpha);	
				for (int run = 0; run < run_max; ++run) {
					long double last_rap = new_rap;
					new_rap = IsoSolver::newRapidity(j, alpha);

					if (finite(new_rap)) break;
					// if it doesn't work, slowly raise the offending rapidity and see if we can get something reasonable.
					rapidity(j, alpha) *= 1.1;
				
					// growing over the machine maximum or otherwise non-finite? give up.
					if (!finite(rapidity(j, alpha))) break; 
				}
				// set the new value reset the old value (don't mess with current rapidities)
				rapidity(j, alpha) = current_rapidity; 
				new_rapidity(j, alpha) = new_rap;
				// if not on hold and not finite, there's nothing we can do.
				if (!finite(new_rapidity(j, alpha))) throw Exception (here, exc_Runaway);
			}
			else { 
				// the value is on hold, don't change.
				new_rapidity(j, alpha) = rapidity(j, alpha);
			}

			// calculate convergence measure
			square_diff += sq(rapidity(j, alpha) - new_rapidity(j, alpha));
			// find deviations & calculate their convergence. if not found, give up.
			if (IsoSolver::getNewDeviances(new_deviance, new_aberration, j, alpha)) {
				for (int a = 0; a < config_.stringLength(j); ++a)
					square_diff += 	  sq( new_deviance(j,alpha,a) - deviance(j,alpha,a) ) 
									+ sq( new_aberration(j,alpha,a) - aberration(j,alpha,a) );
			}
			else {
				throw Exception (here, exc_NonFinite, "no solution for deviance/aberration");
			}
		}
		// if we have NaN square diff, throw exception. *this shouldn't happen.*
		if  (!finite(square_diff)) throw Exception (here, exc_NonFinite); 
				
		// set the new values for the rapidities
		rapidity = new_rapidity;
		deviance = new_deviance;
		aberration = new_aberration;	
		
		// sign off
		convergence = square_diff; 
		++iterations;
		
cerr<<iterations<<' '<<roots()<<endl;
		
	}
	
	
int IsoSolver::sumOfBethe2xJ(const int j, const int alpha) 
{
	// extract sum of Bethe J's from Takahashi I
	int length_j = config_.stringLength(j);
	int sum_2x_bethe_quantum = quantum_number(j, alpha) + chain_length_*(length_j-1);
	for (int k=0; k < config_.numberTypes(); ++k) {
		int length_k = config_.stringLength(k);
		for (int beta=0; beta < config_.numberStringsOfType(k); ++beta) {
			if (j==k && beta==alpha) continue;
			sum_2x_bethe_quantum -= isgn(rapidity(j,alpha) - rapidity(k, beta))*(length_j*length_k - 2*min(length_j, length_k) + (j==k));
		}
	}
	// enforce modulo choice: 2J in (-N..N] (mod 2N)
	while (sum_2x_bethe_quantum > chain_length_) sum_2x_bethe_quantum -= 2*chain_length_;
	while (sum_2x_bethe_quantum <= -chain_length_) sum_2x_bethe_quantum += 2*chain_length_;
	
	return sum_2x_bethe_quantum;
}
		
/** RLH2008 **/				
		
/* find extra real solutions */

bool IsoSolver::findExtraReal(const int j, const int alpha, const int a)
{

	// this should only be called after a collapse. we assume rapidity j, alpha is set to a collapsing solution
	const int max_iter = 1000;
	const double ex_real_precision = 1e-10;
	const double start_distance = 1e-5;

	int jx2;
	if (config_.stringLength(j)==2) {
		// now, we have the sum of two numbers. we also know the two numbers are equal. 
		// so, we know the two numbers, right? right? not right.
		// since we're working modulo N, there's two solutions to that equation. 
		// yes, it pisses me off, too.
		// however, the leap of faith we make is the following: if the numbers were close to zero, 
		// we'd have a regular, christian, tax-paying, real solution. strings are close to the borders.
		// so extra real solutions (which are really just cross-dressing strings) should be found there as well.
		jx2 = sumOfBethe2xJ(j,alpha)/2;
		// I'm also hoping that jx2 is somehow in the interval [-N..N]
		if ( abs(jx2) < chain_length_/2) jx2 -= sgn(jx2)*chain_length_;
		// in which case it's now in the outer halves of that interval. phew!
		// I'm still not quite sure why exactly I need to choose this.
		// however, if we'd pay someone to consider the extremely boring question as to why this is the case, 
		// we might find out.
		// I hate modular arithmetic.
	}
	else {
		// here we need a conjecture for the bethe-J configuration.
		// looking at the sequence of solutions, it looks verry probable that *all* strings consist of sequential quantum numbers.
		return false;
	}
	
cerr<<jx2<<endl;		
	
		
	// first find the trivial root. we only need to check one equation far that.
	// use collapsed rapidity as initial guess
	//TODO: check: we choose the middle rapidity, assuming it has zero aberration. will this work for all cases?
	
	double x_trivial = 0.0;
	double y_trivial = findRoot1(jx2, x_trivial, rapidity(j,alpha));

	if (isNan(y_trivial)) return false;

	// now that we have the trivial root, let's look for the other.
	// here we need to check both equations
	
	// fishing expedition to bracket a second root for positive x
	double x_min = x_trivial + start_distance;
	// find y for this x (stay on root curve for indicator_1), use y_trivial as guess

	double y_min = findRoot1(jx2, x_min, y_trivial);
	if (isNan(y_min)) return false;
	// get the second indicator function
	double f_min = exRealIndicator2(jx2, x_min, y_min);
	
	// outer bracket
	double x_max = x_min;
	double y_max = y_min;
	double f_max;
	
	int iter = 0;

	// fishing expedition for *one* extra root
	do {
		x_max += 1.0;
		// stay on root curve for indicator_1, use last y as guess
		y_max = findRoot1(jx2, x_max, y_max);
		// can't find the curve, drop out
		if (isNan(y_max)) return false;
		// second indicator function
		f_max = exRealIndicator2(jx2, x_max, y_max);

//cout<<iter<<"  ";
//cout<<x_min<<" "<<y_min<<" "<<f_min<<"   ";
//cout<<x_max<<"  "<<y_max<<" "<<f_max<<endl;
		
		if (++iter > max_iter) return false;
	} while ( sgn(f_min) == sgn(f_max) );
//cout<<endl;
	

	double x_try, y_try, f_try;
	// we now have a bracket. zoom in.
	while ( abs(x_max - x_min) + abs(y_max-y_min) >= ex_real_precision ) {
		// interpolate
		//x_try = (abs(f_min)*x_max + abs(f_max)*x_min) / (abs(f_max)+abs(f_min));
		// bisect
		x_try = 0.5*(x_min+x_max);
		// stay on curve. use last y_try as guess
		y_try = findRoot1 (jx2, x_try, y_try); 
		if (isNan(y_try)) return false;
		f_try = exRealIndicator2 (jx2, x_try, y_try);

		if (sgn(f_try)==0) {
			// not bloody likely
			x_min = x_max = x_try;
			y_min = y_max = y_try;
			f_min = f_max = f_try;
		}
		else if (sgn(f_try) == sgn(f_min)) {
			x_min = x_try;	
			y_min = y_try;
			f_min = f_try;
		}
		else {
			x_max = x_try;
			y_max = y_try;
			f_max = f_try;
		}
//cout<<iter<<"  ";
//cout<<x_min<<" "<<y_min<<" "<<f_min<<"   ";
//cout<<x_max<<"  "<<y_max<<" "<<f_max<<endl;
		
		if (++iter > max_iter) return false;
	}

//cout<<indicator1(x_min, y_min)<<" "<< indicator2(x_min, y_min)<<endl;
//cout<<indicator1(x_max, y_max)<<" "<< indicator2(x_max, y_max)<<endl;	
	
	
	// this turns a near pair into two real roots
//FIXME: this only works for two-strings?
	rapidity(j,alpha) = 0.5*(y_min+y_max);
	aberration(j,alpha,a) = 0;
	aberration(j,alpha,a+1) = 0.5*(x_min+x_max);
	deviance(j,alpha,a) = -0.5;
	deviance(j,alpha,a+1) = 0.5;
	
//cout<<roots()<<endl;

	return true;	
}		
	
		
// this is pretty much the same as betheZero, but for extra real solutions.
// (indicator function for a findRoot)
double IsoSolver::exRealIndicator2(const int bethe_2xJ, const double eps, const double lam)
{
 	double scattering_term = 0.0;
 
 	for (int k=0; k < config_.numberTypes(); ++k) 
	for (int beta=0; beta < config_.numberStringsOfType(k); ++beta) 
	for (int b=0; b< config_.stringLength(k); ++b) {
		// we have set the string we're considering separately on hold
		// this is going to be so confusing, aargh.
		// FIXME: we should give j, alpha, a, as parameters and check those instead.
		if (hold_deviation(k,beta,b)) continue;
		
		// we ignore imag parts. they oughta cancel.
		scattering_term += real( atan(lam+eps - root(k, beta,b) ) );
	}	

	return double(chain_length_) * atan(2.0*(lam+eps)) - ( 0.5*PI*bethe_2xJ +  atan(eps) + scattering_term );
}
double IsoSolver::exRealIndicator1(const int bethe_2xJ, const double eps, const double lam)
{
 	double scattering_term = 0.0;

 	for (int k=0; k < config_.numberTypes(); ++k) 
	for (int beta=0; beta < config_.numberStringsOfType(k); ++beta) 
	for (int b=0; b< config_.stringLength(k); ++b) {
		// we have set the string we're considering separately on hold
		// this is going to be so confusing, aargh.
		// FIXME: we should give j, alpha, a, as parameters and check those instead.
		if (hold_deviation(k,beta,b)) continue;
		
		// we ignore imag parts. they oughta cancel.
		scattering_term += real( atan(lam - root(k, beta,b) ) );
	}	
//cerr<<"scat1 "<<scattering_term<<endl;	

	return double(chain_length_) * atan(2.0*lam) - ( 0.5*PI*bethe_2xJ +  atan(-eps) + scattering_term );
}
				
		
double IsoSolver::findRoot1(const int bethe_jx2, const double eps, const double lam_guess)
{
	// find the zero for this eps
	const int max_iter = 5000;
	const double ex_real_precision = 1e-10;
	// not too big as we might skip over the double root
	// TODO: a more sophisticated approach would be more stable here, 
	// e.g. stepping outwards for a number of steps, then looking inwards, etc.
	// but if it's too small then we move too slowly if we don't get it right...
	// i.e we should pick a very small (1e-5?) initial step, then check, then a larger outward step, preferably exponential.
	double fishing_step = 1.0;
	
	
	// use provided lambda as initial guess
	double y_min = lam_guess;
	double y_max = lam_guess;
	double f_min, f_max;

	int iter=0;
	// bracket a root
	do {
		y_max += fishing_step;
		y_min -= fishing_step;
		
		f_min = exRealIndicator1(bethe_jx2, eps, y_min);
		f_max = exRealIndicator1(bethe_jx2, eps, y_max);
		
		if (++iter > max_iter) return NOT_A_NUMBER;
	} while ( sgn(f_min) == sgn(f_max) );


	// we now have a bracket. zoom in.
	while ( abs(y_max - y_min) >= ex_real_precision ) {
			
		double y_try = 0.5*(y_min+y_max);
		double f_try = exRealIndicator1(bethe_jx2, eps, y_try);
		
		if (sgn(f_try)==0) {
			// not bloody likely
			y_min = y_max = y_try;
			f_min = f_max = f_try;
		}
		else if (sgn(f_try) == sgn(f_min)) {
			y_min = y_try;	
			f_min = f_try;
		}
		else {
			y_max = y_try;
			f_max = f_try;
		}
		if (++iter > max_iter) return NOT_A_NUMBER;
	}
	return 0.5*(y_min+y_max);
}








		
		
/** /RLH2008 **/					
		
		
		
		
		
	/** singular ("inadmissible") state solver **/	
	bool IsoSolver::solveSingular(void)
	{
		// sanity check. now, this should be checked before calling; will be harmless&useful here once we store admissibility.
// 		if (admissible()) return false; 	
cerr<<"sing"<<endl;	
		Strip<double> new_rapidity (config_);
		// find the singular pair and put it on hold. more than one: problem.
		int number_singular=0, number_origin=0;
		for (int j=0; j<config_.numberTypes(); ++j) {
			int length_j=config_.stringLength(j);
			if (config_.numberStringsOfType(j) && !quantum_number(j,0)) {
				if (length_j%2) {
					// origin root: we can have only one.
					if (number_origin) return false;
					new_rapidity(j,0) =	0.0;
					hold(j,0) = true;
					hold_deviation(j,0,length_j/2) = true;
					++number_origin;
				}
				else {
					// singular pair: only one.
					if (number_singular) return false;	//throw Exception(here, exc_NoCanDo, "more than one singular pair");
					new_rapidity(j,0) =	0.0;
					hold(j,0) = true;	 
					hold_deviation(j,0,length_j/2-1) = true;	 
					hold_deviation(j,0,length_j/2) = true;	 
					++number_singular;
				}
			}
		}
		
		// iterate the other values
		while (iterations < max_iterations) {
 			iterate();
 			if (abs(convergence) <= precision) return true;
		}
		return false;
	}
	
	
	
	bool IsoSolver::solveSymmetric(void)
	{
// 		if (!oddSymmetric()) return false;	//sanity check
cerr<<" oddsym"<<endl;		
		int number_origin=0;
		int odd_type[2];
		for (int j=0; j<config_.numberTypes(); ++j) {
			int length_j=config_.stringLength(j);
			if (config_.numberStringsOfType(j) && !quantum_number(j,0)) {
				if (length_j%2) {
					// origin root: we can have only two
					if (number_origin>1) return false;
					odd_type[number_origin]=j;
					
					// hold the origin root
					hold(j,0) = true;
					hold_deviation(j,0,length_j/2) = true;
			
					// second odd string: also hold +-i roots.
					if (number_origin) {
						hold_deviation(j,0,length_j/2-1) = true;
						hold_deviation(j,0,length_j/2+1) = true;
					}
					++number_origin;
				}
				else {
					// this state is evil:			
					// singular, odd-symmetric
					// now all bets are off.
					return false;
				}
			}
		}
		int centre_1 = config_.stringLength(odd_type[1])/2;
			
		// both quantum numbers are 0. the Bethe quantum numbers, however, are not.
		// they are \pm 0.5 and \pm (N-1)/2
		double bignm1 = chain_length_-1;
		double oneover_bignm1 = 1.0/bignm1;
		
		// initial values
 		double lambda = 0.01;
 		double delta = 1e-20;
		for (; iterations<max_iterations; ++iterations) {
			complex<double> theta_other = 0.0;
			for (int k=0; k<config_.numberTypes(); ++k)
			for (int beta=0; beta<config_.numberStringsOfType(k);++beta) {
				if (k==odd_type[0] && beta==0) continue;
				if (k==odd_type[1] && beta==0) continue;
				for (int b=0; b<config_.stringLength(k);++b)
					theta_other += atan(lambda - root(k, beta, b)) + atan(I*(1.0+delta) - root(k, beta, b));
			}
				
			
			// dammit, why didn't I think of this before.
			complex<double> lambda_p_imu = tan( 0.5*bignm1*(atan(2.0*lambda) + atan(2.0*I*(1.0+delta))) - 0.5*theta_other);
			double new_lambda = real(lambda_p_imu);
			double new_delta = imag(lambda_p_imu)-1.0;
			
			double sq_diff = sq(new_lambda-lambda) + sq(new_delta-delta);
			
			// set the iteration variables
	  		delta = damping_symmetric*new_delta+(1.0-damping_symmetric)*delta;
  	 		lambda = damping_symmetric*new_lambda+(1.0-damping_symmetric)*lambda;
  	 		
			// and the root data fields
			
			rapidity(odd_type[0],0) = lambda;
			rapidity(odd_type[1],0) = -lambda;
			aberration(odd_type[1],0,centre_1-1) = aberration(odd_type[1],0,centre_1+1) = lambda;
			aberration(odd_type[1],0,centre_1) = 0.0;
			deviance(odd_type[1],0,centre_1-1) = delta;
			deviance(odd_type[1],0,centre_1+1) = -delta;

			
			convergence = 0.0;
			// iterate other roots
			iterate();
			
			
			convergence += sq_diff;
			if (convergence <= precision) return true;
		}
		return false;
	}
			
			
			
			
	

/// TODO: implement the first-order deviation check correctly...
	double IsoSolver::stringDeviation (void) const 
	{

		complex<double> deviation = 0.0;
		double sum_sq = 0;
		for (int j=0; j<config_.numberTypes(); ++j) {
///TODO	: fix this for non-twostrings, for the time being we report them all as deviated...
if (config_.stringLength(j) != 2) return 1;
///
			for (int alpha=0; alpha < config_.numberStringsOfType(j); ++alpha) {
				complex<double> root_j = root(j, alpha, 0);
				deviation = pow(  (root_j + 0.5*I) / (root_j - 0.5*I), -chain_length_);
				// product over other strings
				for (int k=0; k<config_.numberTypes(); ++k)
				for (int beta=0; beta < config_.numberStringsOfType(k); ++beta) {
					if (k==j && beta==alpha) continue;
					for (int b=0; b<config_.stringLength(k);++b){
						complex<double> root_k_b = root(k, beta, b);
						deviation *= (root_j - root_k_b + I) / (root_j - root_k_b - I);
					}
				}
	///TODO	: fix this for non-twostrings
				deviation *= (1.0 + deviance(j, alpha,0));
				deviation -= deviance(j, alpha,0);
	///			
				// norm() gives the *square* of the norm (Bjarne, Bjarne, Bjarne...)
				sum_sq += std::norm(deviation);
			}
		}
		return sqrt(sum_sq);
	}
		
	
	/** gives a rough measure of the amount of deviatedness. for academic purposes. **/
	double IsoSolver::devianceMagnitude(void) const
	{
		double magnitude = 0.0;
		for (int j=0; j < config_.numberTypes(); ++j)
		for (int alpha=0; alpha < config_.numberStringsOfType(j); ++alpha) 
		for (int a=0; a < config_.stringLength(j) ; ++a) 
			magnitude += sq(deviance(j, alpha, a)) + sq(aberration(j, alpha, a));
		return magnitude;
	}

