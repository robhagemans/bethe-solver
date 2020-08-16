/*
            //- zeta(eps1, -del0+del1)
            double lnr_ddel =  (N-1.) * atih(4.+2*del0) - ( 
                              zeta(eps2, 3.+del0) -zeta(eps2,-1.-del0)
                            + zeta(eps1, 4.+del0+del1)  //+zeta(eps1, 2.+del0+del1)-zeta(eps1,-2.-del0-del1)
                            );
            */
            
            
            
        cout<< ( (N-1.)*( xi(2.*eps1,3.+2.*del1) + xi(2.*eps1,-1.-2.*del1) ) - PI*(2.*J+1.) - (
                                   xi(eps1, -del0+del1) + xi(eps1,2.+del0-del1) 
                                 + xi(eps1-eps2, 2.+del1) + xi(eps1-eps2, -del1)
                                 + xi(eps1, 4.+del1+del0) + xi(eps1,-2.-del1-del0) 
                                 + xi(eps1+eps2, 2.+del1) + xi(eps1+eps2, -del1) 
                                 + 2.*atan(2.*eps1) 
                                ))/PI;
        cout<<endl;
        cout<< (N*atan(2.*roots[1]) - atan(roots[1]-roots[0]) - atan(roots[1] - roots[2])  
                                                        - atan(roots[1]-roots[3]) - atan(roots[1] - roots[4]) 
                                                        - atan(roots[1]-roots[5]) - atan(roots[1] - roots[6]) 
                                                        - atan(roots[1]-roots[7]))/PI; 
        cout<<     (N*atan(2.*roots[3]) - atan(roots[3]-roots[0]) - atan(roots[3] - roots[2])  
                                                        - atan(roots[3]-roots[1]) - atan(roots[3] - roots[4]) 
                                                        - atan(roots[3]-roots[5]) - atan(roots[3] - roots[6]) 
                                                        - atan(roots[3]-roots[7]))/PI; 
        
        cout<<endl;
        /*                        
        cout<<  (N-1.) * atan(2.*eps2) - (  
                            + xi(eps2-eps1, 2.+del1) + xi(eps2-eps1, -del1)
                            + xi(eps1+eps2, 2.+del1) + xi(eps1+eps2, -del1)
                            + xi(eps2, 3.+del0)      + xi(eps2, -1.-del0) 
                            );
        cout<<endl;
     
        cout<< N*atan(2.*eps2) - atan(roots[2]-roots[0]) - atan(roots[2] - roots[1])  
                                                        - atan(roots[2]-roots[3]) - atan(roots[2] - roots[4]) 
                                                        - atan(roots[2]-roots[5]) - atan(roots[2] - roots[6]) 
                                                        - atan(roots[2]-roots[7]); 
        */




/*

// determine width of inner pair
		int sum_sign = 0;
		for (int k=0; k < config_.numberTypes(); ++k) 
		for (int beta=0; beta < config_.numberStringsOfType(k); ++beta) {
			if (j==k && alpha==beta) continue;
			int length_k = config_.stringLength(k);
			
			sum_sign += isgn(rapidity(j, alpha)-rapidity(k, beta))*(length_j*length_k - 2*min(length_j, length_k) + (j==k));
		}
		
		int sign_inner = ((quantum_number(j, alpha) + chain_length_*(length_j - 1) - length_j*(config_.numberRoots()+2) - 2 - sum_sign)/2)%2 	? 1 : -1;
*/		
	
/*	
int String::sumJx2(const int j, const int alpha) const 
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
	*/
		


SymmetricSolver::SymmetricSolver(const int N, const std::vector<Roots*> complexes)
    : N_(N), complexes_(complexes)
{
}


/* get all roots in complex form 
*/
vector< complex<double> > SymmetricSolver::getRoots() const
{
 	vector< complex<double> > roots;
    for (int alpha=0; alpha < complexes_.size(); ++alpha) {
        vector< complex<double> > roots_alpha = complexes_[alpha]->getRoots();
        roots.insert(roots.end(), roots_alpha.begin(), roots_alpha.end());
    }
    return roots;
}    

int SymmetricSolver::numberRoots() const
{
 	int number = 0;
    for (int alpha=0; alpha < complexes_.size(); ++alpha) {
        number += complexes_[alpha]->size();
    }
    return number;
}    


/* calculate the quantum numbers (J) as they are in the complex (non-Takahashi) Bethe equations, from the solution 
*/
vector< complex<double> > SymmetricSolver::getDirtyJ() const // (const int N, const vector<Roots*>& strings)
{
    vector< complex<double> > roots = getRoots();
 	vector< complex<double> > qn_J (roots.size());
  	for (int i=0; i < roots.size(); ++i) {
		
        // singular pair needs separate treatment
        // however, we don't need a threshold to check if it's \pm i/2, as these are exact
		if (real(roots[i]) == 0.  &&  abs(imag(roots[i])) == 0.5) {
			int root_sign = isgn(imag(roots[i]));
			if (roots.size()%2) 
				// odd. there's a zero q.n.
				qn_J[i] = root_sign*0.25*(N_ - (N_%4));
			else 
				// even. no zero.
				qn_J[i] = root_sign*0.25*(N_ - (2 - N_%4));
	    }
	    else {	
		    complex<double> lhs = double(N_) * atan(2.*roots[i]);
		    for (int j=0; j<roots.size(); ++j){ 
			    if (j!=i) 
			        lhs -= atan(roots[i] - roots[j]);
		    }
		    qn_J[i] = lhs/PI;
		}	
	}
	return qn_J;
}

vector< double > SymmetricSolver::getRealDirtyJ() const 
{
    vector< complex<double> > roots = getRoots();
 	vector< double > qn_J (roots.size());
  	for (int i=0; i < roots.size(); ++i) {
		
        double lhs_sum = double(N_) * 0.5* (xi(2.*real(roots[i]), 1.+2.*imag(roots[i])) + xi(2.*real(roots[i]), 1.-2.*imag(roots[i])));
	    double lhs_dif = (real(roots[i])==0.) ? double(N_) *rati(2.*imag(roots[i])) : 0.; 
	    
	    for (int j=0; j<roots.size(); ++j) { 
		    if (j!=i) {
		        // add [atan(x) + atan(x*)], if root not real the conjugate will be in the set of roots 
		        lhs_sum -= 0.5*xi(real(roots[i] - roots[j]),1.+imag(roots[i] - roots[j])) ;
		        lhs_sum -= 0.5*xi(real(roots[i] - roots[j]),1.-imag(roots[i] - roots[j])) ;
		        if (real(roots[i]-roots[j])==0.) lhs_dif -= rati(imag(roots[i] - roots[j])) ;
		    }
	    }
	    qn_J[i] = (lhs_sum+lhs_dif)/PI;
		
	}
	return qn_J;
}


vector< double > SymmetricSolver::getImagDirtyJ() const 
{
    vector< complex<double> > roots = getRoots();
 	vector< double > qn_im (roots.size());
  	for (int i=0; i < roots.size(); ++i) {
		
        double lhs_dif = double(N_) * 0.5* (zeta(2.*real(roots[i]), 1.+2.*imag(roots[i])) - zeta(2.*real(roots[i]), 1.-2.*imag(roots[i])));
	    
	    for (int j=0; j<roots.size(); ++j) { 
		    if (j!=i) { 
		        lhs_dif -= 0.5*zeta(real(roots[i] - roots[j]), 1.+imag(roots[i] - roots[j])) ;
		        lhs_dif += 0.5*zeta(real(roots[i] - roots[j]), 1.-imag(roots[i] - roots[j])) ;
		   }
	    }
	    qn_im[i] = lhs_dif;
		
	}
	return qn_im;
}



/* iteration logic 
*/
bool SymmetricSolver::solve(const int max_iter, const double desired_convergence)
{
    for (int i=0; i< complexes_.size(); ++i) {
        complexes_[i]->initiate(N_, complexes_, i);
    }
    
#ifdef showiter   
    cout<<" rts "<< getRoots() <<" ";
    cout<<endl;
#endif

    for (int iter=0; iter<max_iter; ++iter) {

        double convergence = 0.;
        
        for (int i=0; i< complexes_.size(); ++i) {
            complexes_[i]->iterate(N_, complexes_, i, convergence);
        }
        cout.precision(15);

#ifdef showiter    
        cout<< "iter "<< iter <<" ";
        cout<<" conv "<< convergence <<" ";
        cout<<" rts "<< getRoots() <<" ";
        cout<<endl;
#endif

        if (!finite(convergence)) break;        
        if (convergence < desired_convergence) return true;        
    }
    
    // success if converged
    return false;
}



/*
// FIXME: 0.1 is hack, this is only to be used for bethe-takahashi
// xi excluding pi terms, for bethe takahashi
inline double btxi (const double epsilon, const double delta) 
{
    if (delta>0.1) 
    return atan(epsilon/delta);
    else return 0.;
}	
*/
      
        /*
		// scattering with others
		for (int beta=0; beta<all.size(); ++beta) {

		    if (beta==alpha) continue;
		
	        NonsymRoots* roots_beta = dynamic_cast<NonsymRoots*>(all[beta]);
            vector<double> other_x = roots_beta->getRealRoots();
            vector<complex<double> > other_z = roots_beta->getComplexPairs();

    	    for (int k=0; k < other_x.size();++k) {
    	        double x_k = other_x[k];
				scatter += btxi( x-x_k, 1.-y ) + btxi( x-x_k, 1.+y );
		    }
                	
    	    for (int k=0; k < other_z.size();++k) {
    	        double x_k = real(other_z[k]);
    	        double y_k = imag(other_z[k]);
				scatter += btxi( x-x_k, 1.-(y-y_k) )  +  btxi( x-x_k, 1.+(y-y_k) );
				scatter += btxi( x-x_k, 1.-(y+y_k) )  +  btxi( x-x_k, 1.+(y+y_k) );
			}
        }
		*/




/*
		double tan_theta =  tan(  phase - sum_kin  ) ;
		double eps_1 = epsilon_[0];
		double del_1 = delta_[0];
		double del_n2 = delta_[string_length_/2-1];
		
		double b = (eps_1*tan_theta + del_1 - del_n2 + 0.5*string_length_);
		double c = del_n2*(tan_theta*(0.5*string_length_+del_1) - eps_1);
		// keep the sign as it was
		//int sign_sqrt = isgn(2.0*lambda_*tan_theta + b);
		new_lambda = (-b + sign_sqrt*sqrt( sq(b) - 4.0*tan_theta*c ))/(2.0*tan_theta);
*/


/*
    // limit lambda moves - no more than half the distance to nearest rapidity
    if (abs(new_rap-lambda_) > 0.5*smallest_distance) {
        new_rap = lambda_ + 0.499*smallest_distance*sgn(new_rap - lambda_);
        convergence += 1.;
    }
    */


// scale deviations to prevent crossing branch cuts
    // be conservative if lambda moves a lot
    //smallest_distance -= abs(new_rap - lambda_);
    // if (smallest_distance <0.) smallest_distance = 0.;

    for (int a=0; a<string_length_;++a) {
        if (abs(new_epsilon[a]) >= 0.5*smallest_distance) {
        
            
            new_epsilon[a] = 0.49*smallest_distance;
            new_delta[a] = scale*new_delta[a] + (1.-scale)*delta_[a];
            
            convergence+=1.;
        }
    }






/*


    if (has_origin) {
        // 4.L* --> 2* because we're looking at log r^2 rather than log_r, another 2* because we've added eq and conj(eq)
        theta -=  (abs(y.imag())>1.) ? sgn(y.imag()) : 0;
        log_rsq +=  4.L*atih(y.imag()); //log( sq((1.L+y) / (1.L-y)) );
    }
    for (int k=0; k<other_x.size(); ++k) {
        const double x_k = other_x[k];
        // +x // -x
		log_rsq +=  2.L*log( (sq(x_k) + sq(0.5*(2+y.pos)+y.del)) / (sq(x_k) + sq(0.5*(2-y.pos)-y.del)) );
    }
    for (int k=0; k<other_y.size(); ++k) {
        const double del_dif = y.del-other_y[k].del;
        const double del_sum = y.del+other_y[k].del;
        const double pos_dif = 0.5L*(y.pos-other_y[k].pos);
        const double pos_sum = 0.5L*(y.pos+other_y[k].pos);
        // +iy
	    theta -=  (abs(pos_dif+del_dif)>1.) ? sgn(pos_dif+del_dif) : 0;
		log_rsq +=  4.L*atih(pos_dif+del_dif);
        // -iy
	    theta -=  (abs(pos_sum+del_sum)>1.) ? sgn(pos_sum+del_sum) : 0;
		log_rsq +=  4.L*atih(pos_sum+del_sum);
    }

    for (int k=0; k<other_z.size(); ++k) {
        //const long double x_k = other_z[k].real();
        //const long double y_k = imag(other_z[k]);
        // x+iy  -x+iy
        //log_rsq +=  2.L*log( (sq(x_k) + sq(1.L+(y-y_k))) / (sq(x_k) + sq(1.L-(y-y_k))) );
        // -x-iy  x-iy
		//log_rsq +=  2.L*log( (sq(x_k) + sq(1.L+(y+y_k))) / (sq(x_k) + sq(1.L-(y+y_k))) );
        const double re_z = other_z[k].real();
        const int pos_dif = y.pos-other_z[k].pos;
        const int pos_sum = y.pos+other_z[k].pos;
        const double del_dif = y.del - other_z[k].del;
        const double del_sum = y.del + other_z[k].del;
        log_rsq += 2.L*log( (sq(re_z) + sq(0.5*(2+pos_sum)+del_sum)) / (sq(re_z) + sq(0.5*(2-pos_sum)-del_sum)) );
        log_rsq += 2.L*log( (sq(re_z) + sq(0.5*(2+pos_dif)+del_dif)) / (sq(re_z) + sq(0.5*(2-pos_dif)-del_dif)) );
    }

    return;       
*/





vector<complex<double> > CentralString::getCleanerJ (
                            const vector<SymRoots*>& all, const int alpha) const 
{
    // outputs
	vector<complex<double> > result(string_length_);
	
	for (int a=0; a<string_length_; ++a) {
		double log_rsq = 0.0;
	    int theta = 0;
	
	    //const int pos = (string_length_-1) - 2*a;
	    
		// string self scattering 		
		for (int b=0; b < string_length_; ++b) {
			
	        const double y_diff = (delta_[a]-delta_[b]) - (a-b);
	    
			if (b!=a) {    
			    theta -=  (abs(y_diff)>1.) ? sgn(y_diff) : 0;
			    // split position and delta differences for precision in the case where (1-a+b)==0 etc.
    	        log_rsq += log(sq( ((double)(1-a+b) + (delta_[a]-delta_[b])) / ((double)(1+a-b) - (delta_[a]-delta_[b])) ));
            }   
        }		
		
        // scattering with others
        for (int beta=0; beta<all.size(); ++beta) {
            if (beta==alpha) continue;
	          
            const vector<double> other_x = all[beta]->getRealPairs();
            const vector<IVal> other_y = all[beta]->getImagPairs();
            const vector<CVal> other_z = all[beta]->getComplexQuartets();
            Central_scatterOthers( {(string_length_-1) - 2*a, delta_[a]} , 
                        all[beta]->hasOrigin(), other_x, other_y, other_z, log_rsq, theta);
        }   
        
        const double y = 0.5*(string_length_-1) - a + delta_[a];
        // kinetic theta 
		theta += (abs(2.*y)>1.) ? big_n_*sgn(2.*y) : 0;
		log_rsq += (big_n_) * log( sq(y-0.5) / sq(y+0.5));

        // imaginary part should be small, is an error term.
        result[a] = complex<double> (0.5*theta, 0.5*log_rsq);      
    }	
		
    return result;
}

// error is given as average of square error per quantum number
vector<int> CentralString::getCleanJx2 (
                        const vector<SymRoots*>& all, const int alpha, double& sq_error) const
{
    const vector< complex<double> > dirty_j = getCleanerJ(all, alpha);
    vector<int> result(dirty_j.size(), 0);
    sq_error=0.;
    for (int i=0; i<result.size(); ++i) {
        result[i] = round(2.*real(dirty_j[i]));
        sq_error += norm(1.*result[i] - 2.*dirty_j[i])/(double)dirty_j.size();           
    }      
    return result;
}








//    // recalculate Js to account for crossed branch cuts
//        calculateSum2xJ(big_n, all, alpha) ; 
//cerr<<sum_jx2_<<endl;        
        
    /* prevent crossing of branch cuts */
    // if string_length_==1, only place barriers at strings length 3 and over
    // if string_length_==2, only place barriers at strings length 2 and over
    // for longer strings, place barriers at everything.
/*    
    
    // find string extent
    double lambda_left = lambda_;
    double lambda_right = lambda_;
    int furthest_left = string_length_/2;
    int furthest_right = string_length_/2;
    for (int a=0; a<string_length_/2; ++a) {
        if (epsilon_[a] + lambda_ < lambda_left) {
            lambda_left = epsilon_[a] + lambda_;
            furthest_left = a;
        }
        if (epsilon_[a] + lambda_ > lambda_right) {
            lambda_right = epsilon_[a] + lambda_;
            furthest_right = a;
        }
    }
    
    // find nearest other rapidity
    double smallest_distance_right = big_n_;
    double smallest_distance_left = big_n_;
    
    for (int beta=0; beta<all.size(); ++beta) {
        if (beta==alpha) continue;
        vector<CVal> other_lambda = (all[beta])->getRoots();
        
        if (string_length_==1 && all[beta]->stringLength()<3) continue;
        if (string_length_==2 && all[beta]->stringLength()==1) continue;
        
        
        for (int k=0; k < other_lambda.size();++k) {
            const double x_k = other_lambda[k].real();
            
            if (x_k<lambda_left && lambda_left-x_k<smallest_distance_left) 
                smallest_distance_left = lambda_left-x_k;
            if (x_k>lambda_right && x_k-lambda_right<smallest_distance_right) 
                smallest_distance_right = x_k-lambda_right;
        }    
    }

//#ifdef full_limits        
   
    //if (smallest_distance==0.) return false;
    bool limit_left = false;
    bool limit_right = false;
    // moving right
    if (new_lambda_-lambda_ > 0.5*smallest_distance_right) {
        new_lambda_ = lambda_ + 0.499*smallest_distance_right;
        convergence += 1.;
        limit_right = true;
        result = IterResult::iter_constrained;
    }
    
    // moving left
    if (lambda_-new_lambda_ > 0.5*smallest_distance_left) {
        new_lambda_ = lambda_ - 0.499*smallest_distance_left;
        
        convergence += 1.;
        limit_left = true;
        result = IterResult::iter_constrained;

    }
//#endif


*/




/*
//#ifdef full_limits
    

    if (limit_left) {
        for (int a=0; a<string_length_;++a) {
            if (new_epsilon[a] < epsilon_[furthest_left]) {
                new_epsilon[a] = epsilon_[furthest_left];
                convergence += 1.;
                result = IterResult::iter_constrained;
            }
        }
    }
    else if (limit_right) {
        for (int a=0; a<string_length_;++a) {
            if (new_epsilon[a] > epsilon_[furthest_right]) {
                new_epsilon[a] = epsilon_[furthest_right];
                convergence += 1.;
                result = IterResult::iter_constrained;
            }
        }
    }
    else 
//#endif

    {
        for (int a=0; a<string_length_;++a) {
            // moving left
            if (new_epsilon[a] < 0. && new_lambda_ + new_epsilon[a] - lambda_ - epsilon_[a] <= -0.5*smallest_distance_left) {
                new_epsilon[a] = epsilon_[a] - 0.49*smallest_distance_left + (lambda_ - new_lambda_);
                convergence +=1.;
                result = IterResult::iter_constrained;

            }
            // moving right
            else if (new_epsilon[a] > 0. && new_lambda_ + new_epsilon[a] - lambda_ - epsilon_[a] >= 0.5*smallest_distance_right) {
                new_epsilon[a] = epsilon_[a] + 0.49*smallest_distance_right + (lambda_ - new_lambda_);
                convergence +=1.;
                result = IterResult::iter_constrained;

            }
        }
        
    }            
*/    
