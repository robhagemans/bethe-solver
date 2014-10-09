
#include "roots.h"

using std::vector;
using std::complex;
using std::cout;
using std::cerr;
using std::endl;
using std::abs;



// values for vector of CVal
std::vector< std::complex<double> > getValues(const std::vector<CVal>& roots)
{
 	std::vector< std::complex<double> > root_vals (roots.size());
    for (int i=0; i < roots.size(); ++i) {
        root_vals[i] = roots[i].val();
    }
    return root_vals;
}    



// bethe quantum numbers
vector< complex<double> > getJs(const int big_n, const vector<CVal>& roots) 
{
 	vector< complex<double> > qn_j (roots.size());
  	for (int i=0; i < roots.size(); ++i) {
		
        // singular pair needs separate treatment
        if ( roots[i].real() == 0.  &&  abs(roots[i].pos)==1 ) 
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
	    {	
		    //double sum_theta = 0.5*big_n* re_atan_x2(mul(roots[i], 2));
		    HPIVal sum_theta = mul( re_atan(mul(roots[i], 2)), 2*big_n);
		    double sum_log_rsq = 0.25*big_n* im_atan_x4(mul(roots[i], 2));

		    for (int j=0; j<roots.size(); ++j){ 
			    if (j==i) continue;
		     
		        //sum_theta.sub(re_atan_x2( dif(roots[i], roots[j]) ));
		        sum_theta.sub(mul(re_atan( dif(roots[i], roots[j]) ),2));
		        
		        sum_log_rsq -= 0.25*im_atan_x4(dif(roots[i], roots[j]));     
		    }
		    //qn_j[i] = (1./PI)*(sum_theta.val()*0.5 + I*sum_log_rsq);
            qn_j[i] = M_1_PI*(sum_theta.val()*0.5 + I*sum_log_rsq);
        
        }
    }	
	return qn_j;
}


vector<int> getJx2 (const int big_n, const vector<CVal>& roots, double* error_per_number, bool* parity_error) 
{
    const vector< complex<double> > dirty_j = getJs(big_n, roots);
    int n_roots = dirty_j.size();
    bool jx2_even = (big_n + n_roots)%2;
    vector<int> jx2(n_roots, 0);
    
    if (error_per_number) 
        *error_per_number = 0.;
    if (parity_error) 
        *parity_error = false;
    
    for (int i=0; i<n_roots; ++i) {
        jx2[i] = round(2.*real(dirty_j[i]));
        
        if (parity_error && jx2[i]%2 == jx2_even) 
            *parity_error = true; 
        if (error_per_number)
            *error_per_number += norm((double)jx2[i] - 2.*dirty_j[i]);           
    }      
    
    if (error_per_number) 
        *error_per_number = sqrt(*error_per_number/n_roots);
    
    return jx2;
}



// bethe quantum numbers, sanity check
// some cancellations will not work out accurate so you'll end up with imaginary parts if deltas are small
vector< complex<double> > getDirtyJs(const int big_n, const vector<CVal>& roots) 
{
    vector< complex<double> > qn_j (roots.size());
  	for (int i=0; i < roots.size(); ++i) {
        complex<double> lhs = double(big_n) * std::atan(2.*roots[i].val());
	    for (int j=0; j<roots.size(); ++j){ 
		    if (j!=i)  
		        lhs -= std::atan(complex<double> (
		            (roots[i].real()-roots[j].real()), 0.5*(roots[i].pos-roots[j].pos) + (roots[i].del-roots[j].del) ));
	    }
	    qn_j[i] = (lhs/PI);
    }	
	return qn_j;
}





/** symmetric scattering **/


HPIVal sum_re_atan(const double x,  
        const bool has_origin, const vector<double>& other_x, 
        const vector<IVal>& other_y, const vector<CVal>& other_z) 
{
    HPIVal theta = { 0, 0 };
    if (has_origin) {
        theta.add( re_atan(x) );
    }
    for (int k=0; k< other_x.size(); ++k) { 
        theta.add( sum(re_atan(x-other_x[k]), re_atan(x+other_x[k])) );
        // WARNING - low precision for string differences close to zero 
        // as we're not differentiating lambda and epsilon
        // RVal could help fix this but we'd also need to introduce nonzero epsilons on the real axis
        // perhaps for N=12 , -8, -4 state?
    }    
    for (int k=0; k< other_y.size(); ++k) {
//        theta += 0.5*sum(re_atan_x2(dif(x, other_y[k])), re_atan_x2(sum(x, other_y[k]))).val();
        theta.add( sum(re_atan(dif(x, other_y[k])), re_atan(sum(x, other_y[k]))) );

    }
    for (int k=0; k< other_z.size(); ++k) {
//        theta += sum(re_atan_x2(dif(x, other_z[k])), re_atan_x2(sum(x, other_z[k]))).val();
        theta.add(mul( sum(re_atan(dif(x, other_z[k])), re_atan(sum(x, other_z[k]))), 2));
    }
    return theta;
}


double sum_im_atan(const IVal y,  
        const bool has_origin, const vector<double>& other_x, 
        const vector<IVal>& other_y, const vector<CVal>& other_z, 
        int& itheta)  
{
    double log_r = 0.;
    if (has_origin) {
        log_r += 0.25*im_atan_x4(y);
        //itheta += re_atan_x4_div_pi(y);
        itheta += 2*re_atan_x2_div_pi(y);
    
    }
    for (int k=0; k< other_x.size(); ++k) {
        log_r += 0.5*im_atan_x4(dif(y, other_x[k])); 
    }
    for (int k=0; k< other_y.size(); ++k) {  
        log_r += 0.25*(im_atan_x4(dif(y, other_y[k])) + im_atan_x4(sum(y, other_y[k])));
        //itheta += re_atan_x4_div_pi(dif(y, other_y[k])) + re_atan_x4_div_pi(sum(y, other_y[k]));
        itheta += 2*(re_atan_x2_div_pi(dif(y, other_y[k])) + re_atan_x2_div_pi(sum(y, other_y[k])));
    }    
    for (int k=0; k< other_z.size(); ++k) {
        log_r += 0.5*(im_atan_x4(sum(y, other_z[k])) + im_atan_x4(dif(y, other_z[k])));     
    }    
    return log_r;
}

//Sym_thetaOthers
HPIVal sum_re_atan (const CVal& z, 
                            const bool has_origin, const vector<double>& other_x, 
                            const vector<IVal>& other_y, const vector<CVal>& other_z) 
{
 	HPIVal theta = { 0, 0. };
	if (has_origin) {
        theta.add( re_atan (z) );
    }
    for (int k=0; k<other_x.size(); ++k) {
        theta.add(re_atan( dif(z, other_x[k]) ));
        theta.add(re_atan( sum(z, other_x[k]) ));
    }
    for (int k=0; k<other_y.size(); ++k) {
        theta.add(re_atan( dif(z, other_y[k]) ));     
        theta.add(re_atan( sum(z, other_y[k]) ));
    }
    for (int k=0; k<other_z.size(); ++k) {
        theta.add(re_atan( dif(z, other_z[k]) ));
        theta.add(re_atan( sum(z, other_z[k]) ));
        theta.add(re_atan( dif(z, conj(other_z[k])) ));
        theta.add(re_atan( sum(z, conj(other_z[k])) ));
    }
    return theta;
}

//Sym_logRSqOthers
double sum_im_atan_x4 (const CVal& z, 
                            const bool has_origin, const vector<double>& other_x, 
                            const vector<IVal>& other_y, const vector<CVal>& other_z) 
{
	double log_rsq = 0.;
    if (has_origin) {
        log_rsq += im_atan_x4(z);    
    }
    for (int k=0; k<other_x.size(); ++k) {
        log_rsq += im_atan_x4(dif(z, other_x[k])) + im_atan_x4(sum(z, other_x[k]));
    }
    for (int k=0; k<other_y.size(); ++k) {
        log_rsq += im_atan_x4(dif(z, other_y[k])) + im_atan_x4(sum(z, other_y[k]));
    }
    for (int k=0; k<other_z.size(); ++k) {
        log_rsq += im_atan_x4(dif(z, other_z[k])) + im_atan_x4(sum(z, other_z[k]));
        log_rsq += im_atan_x4(dif(z, conj(other_z[k]))) + im_atan_x4(sum(z, conj(other_z[k])));
    }
    return log_rsq;
}


/** nonsymmetric scattering **/


// real root
HPIVal sum_re_atan (const double x, const vector<double>& other_x, const vector<CVal>& other_z) 
{
    HPIVal theta = { 0, 0 };  
    for (int k=0; k<other_x.size(); ++k) {
        theta.add(re_atan( x-other_x[k] ));
    }
    for (int k=0; k<other_z.size(); ++k) {
        theta.add(mul(re_atan(dif( x, other_z[k] )) ,2));
    }
    return theta;    
}


// complex root, real part
HPIVal sum_re_atan (const CVal& z, const vector<double>& other_x, const vector<CVal>& other_z) 
{
    HPIVal theta = { 0, 0. };  
    for (int k=0; k<other_x.size(); ++k) {
        theta.add(re_atan( dif(z, other_x[k]) ));
    }
    for (int k=0; k<other_z.size(); ++k) {
        theta.add(re_atan( dif(z, other_z[k]) ));
		theta.add(re_atan( dif(z, conj(other_z[k])) ));
    }
    return theta;    
}

// complex root, imag part
double sum_im_atan_x4 (const CVal& z, const vector<double>& other_x, const vector<CVal>& other_z) 
{
    double log_rsq = 0.;
    for (int k=0; k<other_x.size(); ++k) {
        log_rsq += im_atan_x4( dif(z, other_x[k]) );
    }
    for (int k=0; k<other_z.size(); ++k) {
        log_rsq += im_atan_x4( dif(z, other_z[k]) );
		log_rsq += im_atan_x4( dif(z, conj(other_z[k])) );
    }
    return log_rsq;
}



        
