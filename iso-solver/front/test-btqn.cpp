
#include "roots.h"
#include "nonsymmetric-string.h"
#include "solver.h"


#include "roots.h"
#include "symmetric-string.h"
#include "solver.h"

using namespace std;



/** for use in Bethe-Takahashi approximation **/

inline HPIVal xi_plus_BT (const CVal& z)
{   return neg(im_log_i(add_n_i(z, 1)));   }

inline HPIVal xi_minus_BT (const CVal& z)
{   return im_log_i(add_n_i(z, -1));   }

inline HPIVal re_atan_BT (const CVal& z)
{   return div(sum(xi_plus_BT(z), xi_minus_BT(z)),2);  }    


//inline CVal pos_div_plus(const CVal& z)
//{   return { z.lam/(1-z.imag()), z.eps)/(1-z.imag()), 0, 0 };  } 

inline HPIVal re_atan_plus (const CVal& z)
{
    //HPIVal reat = (z.pos==2 ? re_atan( {z.lam/-z.del,z.eps/-z.del,0,0} ) : re_atan(pos_div_plus(z)));
//std::cerr<<z.val()<<" "<<z.pos<<" "<<pos_div_plus(z).val()<<" at "<<reat.val()<<endl;
    if (z.pos==2) return { 0, re_atan(z.real()/(0.5*(2-z.pos)-z.del)).eps }; // thow away pi/2 if it's there, absorb into B_T qn.
    return re_atan(z.real()/(0.5*(2-z.pos)-z.del));
}
//return (z.pos==2 ? re_atan( {z.lam/-z.del,z.eps/-z.del,0,0} ) : re_atan(pos_div_plus(z))); }    


// bethe-takahashi quantum numbers
// roots must be organised in vector of vectors with each inner vector corresponding to a string
std::vector<std::complex<double> > getIs_BT(const int big_n, const std::vector<std::vector<CVal> >& roots);


// bethe quantum numbers
vector< complex<double> > getJs_BT_OLD(const int big_n, const vector<CVal>& roots) 
{
    //vector< complex<double> > roots = getRoots();
 	vector< complex<double> > qn_j (roots.size());
  	for (int i=0; i < roots.size(); ++i) {
		
        // singular pair needs separate treatment
/*        if ( roots[i].real() == 0.  &&  abs(roots[i].pos)==1 ) 
		{
		   
		    const int root_sign = isgn(roots[i].pos);
			if (roots.size()%2) 
				// odd. there's a zero q.n.
				qn_j[i] = root_sign*0.25*(big_n - (big_n%4));
			else 
				// even. no zero.
				qn_j[i] = root_sign*0.25*(big_n - (2 - big_n%4));
	    }
	    else 
*/	    {	
		    HPIVal sum_theta = mul( re_atan_BT(mul(roots[i], 2)), big_n);
		    double sum_log_rsq = big_n*im_atan_x4(mul(roots[i], 2));
		    
		    for (int j=0; j<roots.size(); ++j){ 
			    if (j==i) continue;
		     
		        sum_theta.sub(re_atan_BT( dif(roots[i], roots[j]) ));
		        
		        HPIVal s = neg( re_atan_BT( dif(roots[i], roots[j]) ) );
    	        sum_log_rsq -= im_atan_x4(dif(roots[i], roots[j]));     
		    }

		    qn_j[i] = M_1_PI*(sum_theta.val() + I*0.25*sum_log_rsq);
        }
    }	
	return qn_j;
}


// bethe-takahashi quantum numbers
// roots must be organised in vector of vectors with each inner vector corresponding to a string
vector<complex<double> > getIs_BT_OLD(const int big_n, const vector<vector<CVal> >& roots) 
{
    vector< complex<double> > qn_j (roots.size());
  	for (int alpha=0; alpha < roots.size(); ++alpha) {
  	    const int n_alpha = roots[alpha].size();
        HPIVal sum_theta = {0,0};
        double sum_log_rsq = 0;
        for (int a=0; a<roots[alpha].size(); ++a) {
        		    
		    sum_theta.add( mul( re_atan_BT(mul(roots[alpha][a], 2)), big_n) );
		    sum_theta.sub_half_pi( big_n*(n_alpha-1)*isgn(roots[alpha][a].real()) );
		    sum_log_rsq += big_n*im_atan_x4(mul(roots[alpha][a], 2));
		    
		    for (int beta=0; beta<roots.size(); ++beta) {
		        const int n_beta = roots[beta].size();
		        for (int b=0; b<roots[beta].size(); ++b)
		        { 
			        if (beta==alpha && a==b) continue; 
			        
		            sum_theta.sub( re_atan_BT(dif(roots[alpha][a], roots[beta][b])) );
		            sum_log_rsq -= im_atan_x4(dif(roots[alpha][a], roots[beta][b]));     
		        }
		        // subtract pi/2 for each atan term in Bethe-Takahashi equations
		        sum_theta.add_half_pi( isgn(dif(roots[alpha][a],roots[beta][0]).real())*( -1 + (n_alpha==n_beta) + 2*std::min(n_alpha, n_beta)) );
		    }
        }
        qn_j[alpha] = M_1_PI*(sum_theta.val() + I*0.25*sum_log_rsq);
    }	
	return qn_j;
}


// bethe-takahashi quantum numbers
// roots must be organised in vector of vectors with each inner vector corresponding to a string
vector<complex<double> > getIs_BT(const int big_n, const vector<vector<CVal> >& roots) 
{
    vector< complex<double> > qn_j (roots.size());
  	for (int alpha=0; alpha < roots.size(); ++alpha) {
  	    const int n_alpha = roots[alpha].size();
        HPIVal sum_theta = {0,0};
        double sum_log_rsq = 0;
        for (int a=0; a<roots[alpha].size(); ++a) {
cerr<<endl<<"K "<<alpha<<" "<<a<<" ";
cerr<<mul(roots[alpha][a], 2).val()<<" "<<re_atan_plus(mul(roots[alpha][a], 2)).val()<<" ";
cerr<<re_atan(mul(roots[alpha][a], 2)).val()<<" ";
cerr<<std::atan(mul(roots[alpha][a], 2).val())<<" ";
cerr<<endl;  
		    sum_theta.add( mul( re_atan_plus(mul(roots[alpha][a], 2)), big_n) );
		    sum_log_rsq += big_n*im_atan_x4(mul(roots[alpha][a], 2));
		    
		    for (int beta=0; beta<roots.size(); ++beta) {
		        const int n_beta = roots[beta].size();
		        for (int b=0; b<roots[beta].size(); ++b)
		        { 
			        if (beta==alpha && a==b) continue; 
cerr<<"S "<<alpha<<" "<<a<<"   "<<beta<<" "<<b<<endl;			       
cerr<<dif(roots[alpha][a], roots[beta][b]).val()<<" "<<re_atan_plus(dif(roots[alpha][a], roots[beta][b])).val()<<" ";
cerr<<re_atan(dif(roots[alpha][a], roots[beta][b])).val()<<" ";
cerr<<std::atan(dif(roots[alpha][a], roots[beta][b]).val())<<" ";
cerr<<endl;
			        
		            sum_theta.sub( re_atan_plus(dif(roots[alpha][a], roots[beta][b])) );
		            sum_log_rsq -= im_atan_x4(dif(roots[alpha][a], roots[beta][b]));     
		        }
		        
		    }
        }
        qn_j[alpha] = M_1_PI*(sum_theta.val() + I*0.25*sum_log_rsq);
    }	
	return qn_j;
}


void showIteration (StateSolver* solver)
{
    cout<< "iter "<< solver->iterations() <<" ";
    cout<<" conv "<< solver->convergence() <<" ";
    cout<<" rts "<< getValues(solver->getRoots()) <<" ";
    cout<<endl;
};



int main()
{
    double desired_convergence = 1e-20;
    
    int N;
    StateSolver* all_roots = 0;
    
    
    int max_iter = 200;
    
    //for(int test=1; test<=100;++test) 
    {
           
            /*
             state [1 1 1]#19
 2I [[-4] [-2] [+0]] 
-3
12
-3
 iter 49
 conv 6.56439e-21
 roots [(-0.303681,0) (-0.641278,0.50335) (-0.641278,-0.50335) (0.529579,0.996605) (0.527081,0) (0.529579,-0.996605)]
 jx2 [-3 -7 -5 7 5 9]
 err 1.57954e-10

            */
            
        N=12;
        all_roots = new Bunch<NonsymRoots> ({ new NonsymString(N, 1, -4) ,new NonsymString(N, 2, -2), new NonsymString(N, 3, 0) }, &showIteration);


        if (!all_roots) {
            cout<<" no test case with that number"<<endl;
            return 1;
        }
        
        IterResult solved = all_roots->solve(max_iter, desired_convergence);
        vector<CVal> roots = all_roots->getRoots();
        delete all_roots;
        
        
        vector<CVal> roots1 = roots;//all_roots->getRoots();
            vector< vector<CVal> > roots_string = {
                vector<CVal> { roots1[0] },
                vector<CVal> { roots1[1], roots1[2] },
                vector<CVal> { roots1[3], roots1[4], roots1[5] }
                };
            cout<<"I_BT? "<<getIs_BT(N, roots_string)<<endl;            
            
        
        cout.precision(15);
//        cout<<endl<<" test "<<test<<endl;
        cout<<" N "<<N<<endl; 
        cout<<" roots "<<getValues(roots)<<endl;
        cout.precision(4);
        cout<<endl<<" dirty "<< getDirtyJs(N, roots) <<endl;
        cout<<endl<<" clean "<< getJs(N, roots) <<endl;
//        cout<<endl<<" BT "<< getJs_BT(N, roots) <<endl;

        double betherr = 0;
        bool parityerr = false;
        cout<<endl<<" jx2 "<< getJx2(N, roots, &betherr, &parityerr) <<endl;
        cout<<" err "<< betherr<<endl;
        cout<<endl;
        if (solved != IterResult::iter_ok) 
            cout<<"\033[1;31m "<<message(solved)<<" \033[0m"<<endl; 
        if (parityerr) 
            cout<<"\033[1;31m parity error \033[0m"<<endl; 
        if (betherr>1e-5)
            cout<<"\033[1;31m large quantum number error \033[0m"<<endl;
    }
    return 0;    
}




