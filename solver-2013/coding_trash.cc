bool XXXDeviatedState::findExtraReal(const int j, const int alpha, const int a)
{

	// this should only be called after a collapse. we assume rapidity j, alpha is set to a collapsing solution
	const int max_iter = 1000;
	const double ex_real_precision = 1e-10;
	const double start_distance = 0.2; // 1e-5

	
	// extract sum of Bethe J's from Takahashi I
	int length_j = p_chain->stringLength(j);
	int sum_2x_bethe_quantum = quantum_number(j, alpha) + p_chain->length()*(length_j-1);
	for (int k=0; k < p_base->numberTypes(); ++k) {
		int length_k = p_chain->stringLength(k);
		for (int beta=0; beta < p_base->numberStringsOfType(k); ++beta) {
			if (j==k && beta==alpha) continue;
			sum_2x_bethe_quantum -= isgn(rapidity(j,alpha) - rapidity(k, beta))*(length_j*length_k - 2*min(length_j, length_k) + (j==k));
		}
	}
	
	// initial brackets, to be narrowed down by nearest branch cuts
	double barrier_min = -1000;
	double barrier_max = 1000;
	
	int jx2;
	if (length_j==2) {
		// now, we have the sum of two numbers. we also know the two numbers are equal. 
		// since we're working modulo N, there's two solutions to that equation: 2J = sum (mod N)  ==>  J = sum/2 (mod N/2) 
        // we expect to find string-like solutions near the border of the accepted interval.
		jx2 = sum_2x_bethe_quantum/2;
		// I'm also hoping that jx2 is somehow in the interval [-N..N]
		if ( abs(jx2) < p_chain->length()/2) jx2 -= sgn(jx2)*p_chain->length();
		// in which case it's now in the outer halves of that interval. phew!
		// I'm still not quite sure why exactly I need to choose this.
	}
	else {
	    // real part of current root. to establish initial brackets to avoid crossing branch cuts.
		double real_root = real(root(j,alpha,a));
		
		
		// here we need a conjecture for the bethe-J configuration.
	    // we assume the quantum numbers for the collapsed solution are correct
	    vector<int> all_jx2 = calculateBethe2I();

	    // get the quantum number corresponding to (j,alpha,a) out of this vector:
	    // this should be equivalent to that for (j,alpha,a+1) as these are a collapsed pair.
	    // also establish barrier brackets in this loop at the location of the nearest string rapidities.
	    for (int index=0, k=0; k < p_base->numberTypes(); ++k)  
	    for (int beta=0; beta < p_base->numberStringsOfType(k); ++beta) 
		for (int b=0; b < p_chain->stringLength(k); ++b, ++index) {

		    if ((k==j) && (beta==alpha) && (b==a))  {
		        jx2 = all_jx2[index];
		    }
		    else {
		    
/** SHIFT **  if (b!=a+1) {
                cerr<<"dev "<<deviance(k,beta,b)<<endl;
                  deviance(k,beta,b) = -0.35*deviance(k,beta,b);
               }
/** SHIFT  **/	
	
	
		        complex<double> root_k = root(k,beta,b);
		        if (imag(root_k) > 1e-8) {
		            double real_root_k = real(root_k); 
		        
		            if ((real_root_k < real_root) && (barrier_min < real_root_k))  barrier_min = real_root_k;
		            if ((real_root_k > real_root) && (barrier_max > real_root_k))  barrier_max = real_root_k;
		        }    
		    }
		}  
	}
cerr<<"roots: "<<roots()<<endl;
cerr<<"jx2: "<< jx2<<endl;		
cerr<<"barriers: "<<barrier_min<<" "<<barrier_max<<endl;	
		
		
   cerr<<endl<<endl;
   double x_from = -50;
   double x_to = 50;

   double y_from = barrier_min-50 ;
   double y_to = barrier_min +50;
   
   int x_steps = 65;
   int y_steps = 120;
   
   cerr << x_from << " to " << x_to <<endl;
   for(int iy=0; iy < y_steps; ++iy) {
       
       double y_try = y_from +  iy * (y_to-y_from)/(1.0*y_steps); // position
       double f_try = 0.0;
       
       for(int ix=0; ix < x_steps; ++ix) {
           double x_try = x_from + ix * (x_to - x_from)/(1.0*x_steps); // distance
           
           //double last_f = f_try;
           f_try =  exRealIndicator1(jx2, x_try, y_try) * exRealIndicator2(jx2, x_try, y_try);
           double phase = 2.0* f_try / PI;
           //if (last_f*f_try<0) cerr<<x_try;
           if (f_try>0) cerr <<"+"; else cerr<<"-";
           cerr<< char(48+ int(abs(floor(phase))));

     //      if (exRealIndicator1(jx2, x_try, y_try)>0) cerr <<"+"; else cerr<<"-";
     //      if (exRealIndicator2(jx2, x_try, y_try)>0) cerr <<"1"; else cerr<<"0";
           
       }
       cerr<<" "<<y_try<<endl;
   }

   return false;
   

		
		
	// first find the trivial root. we only need to check one equation far that.
	// use collapsed rapidity as initial guess
	
	for (int it = 0; it< 10; ++it) {
	    double x_trivial = 0.0;
	    double y_trivial = findRoot1(jx2, x_trivial, rapidity(j,alpha), barrier_min);

    cerr<< "y trivial "<< y_trivial<<" f_trivial "<<exRealIndicator2(jx2, x_trivial, y_trivial)<<endl;
	    if (isNan(y_trivial)) return false;
        

		
	    // now that we have the trivial root, let's look for the other.
	    // here we need to check both equations

	    // fishing expedition to bracket a second root for positive x
	    double x_min = x_trivial + start_distance; //// pos x 

	    // find y for this x (stay on root curve for indicator_1), use y_trivial as guess
	    double y_min = findRoot1(jx2, x_min, y_trivial, barrier_min);
	    if (isNan(y_min)) return false;
	
	
	
	
	    // this turns a near pair into two real roots
	    rapidity(j,alpha) = y_min;
	    aberration(j,alpha,a) = 0;
	    aberration(j,alpha,a+1) = x_min;
	    deviance(j,alpha,a) = -0.5;
	    deviance(j,alpha,a+1) = 0.5;
	    iterate();
      
    cerr<<"roots: "<<roots()<<endl;
    }
    
		
    
return false;
	
/*	
	// get the second indicator function
	double f_min = exRealIndicator2(jx2, x_min, y_min);
	
	// outer bracket
	double x_max = x_min;
	double y_max = y_min;
	double f_max;
	
	int iter = 0;

	// fishing expedition for *one* extra root
	do {
		//x_max += 0.1;
        x_max += 0.1;

		// stay on root curve for indicator_1, use last y as guess
		y_max = findRoot1(jx2, x_max, y_max, barrier_min);
		// can't find the curve, drop out
		if (isNan(y_max)) return false;
		// second indicator function
		f_max = exRealIndicator2(jx2, x_max, y_max);

cerr<<iter<<"  ";
cerr<<x_min<<" "<<y_min<<" "<<f_min<<"   ";
cerr<<x_max<<"  "<<y_max<<" "<<f_max<<endl;
		
		if (++iter > max_iter) return false;
	} while ( sgn(f_min) == sgn(f_max) );
cerr<<endl;
	

	double x_try, y_try, f_try;
	// we now have a bracket. zoom in.
	while ( abs(x_max - x_min) + abs(y_max-y_min) >= ex_real_precision ) {
		// interpolate
		//x_try = (abs(f_min)*x_max + abs(f_max)*x_min) / (abs(f_max)+abs(f_min));
		// bisect
		x_try = 0.5*(x_min+x_max);
		// stay on curve. use last y_try as guess
		y_try = findRoot1 (jx2, x_try, y_try, barrier_min); 
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
cout<<"ZOOM: "<<iter<<"  ";
cout<<x_min<<" "<<y_min<<" "<<f_min<<"   ";
cout<<x_max<<"  "<<y_max<<" "<<f_max<<endl;
		
		if (++iter > max_iter) return false;
	}
*/	
	
	
/*	
	// this turns a near pair into two real roots
	rapidity(j,alpha) = 0.5*(y_min+y_max);
	aberration(j,alpha,a) = 0;
	aberration(j,alpha,a+1) = 0.5*(x_min+x_max);
	deviance(j,alpha,a) = -0.5;
	deviance(j,alpha,a+1) = 0.5;
*/

	
	return true;	
}		

