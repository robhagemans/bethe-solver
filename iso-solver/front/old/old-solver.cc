


class Solver {
public:
    Solver(const int big_n, const std::vector<Roots*> complexes);
    Solver(const int big_n, const std::vector<int>& base, const std::vector< std::vector<int> >& ix2);

    std::vector< std::complex<double> > getRoots() const;
    
    std::vector< std::complex<double> > getDirtyJ() const;
    int numberRoots() const;

    template<class ReportFunc = NoFunc>
    bool solve(const int max_iter, const double desired_convergence, ReportFunc report = NoFunc());
     
private:
    // construction succeeded
    bool good_;
    // chain length
    int big_n_;

    // all root complexes    
    std::vector<Roots*> complexes_;


public:
    int iter_;
    double conv_;
      
};


/* iteration logic 
*/

template<class ReportFunc>
bool Solver::solve(const int max_iter, const double desired_convergence, ReportFunc report)
{
    if (!good_) return false;
    
    for (int i=0; i< complexes_.size(); ++i) {
        complexes_[i]->initiate(big_n_, complexes_, i);
    }
    report(this);
    
    for (; iter_<max_iter; ++iter_) {

        conv_ = 0.;
        
        // calculate next value
        for (int i=0; i< complexes_.size(); ++i) {
            if (!(complexes_[i]->iterate(big_n_, complexes_, i, conv_, conv_))) return false;
        }
        // replace old value with new
        for (int i=0; i< complexes_.size(); ++i) {
            complexes_[i]->refresh();
        }
        
        report(this);
        
        if (!finite(conv_)) break;        
        if (conv_ < desired_convergence) return true;        
    }
    
    // success if converged
    return false;
}




Solver::Solver(const int big_n, const std::vector<Roots*> complexes)
    : good_(true), big_n_(big_n), complexes_(complexes), iter_(0), conv_(1.)
{
}


// string configuration constructor
Solver::Solver(const int big_n, const vector<int>& base, const vector< vector<int> >& ix2)
    : good_(true), big_n_(big_n), complexes_(), iter_(0), conv_(1.)
{

    if (!isSymmetric(ix2)) {
        for(int j=0; j<ix2.size(); ++j) 
        for(int alpha=0; alpha<ix2[j].size(); ++alpha) 
            complexes_.push_back(  new NonsymString( length(j), ix2[j][alpha] )  );
        
        return; 
    }
    
    vector<int> pos_real;
    vector<int> odd_zeros;
    vector<int> even_zeros;
    
    for(int alpha=0; alpha<ix2.at(0).size(); ++alpha) {
        if (ix2[0][alpha] == 0) 
            odd_zeros.push_back(1);
            //result.push_back( new OriginRoot() );
        
        if (ix2[0][alpha] > 0) 
            pos_real.push_back(guessSum2xJ(big_n, ix2, 0, alpha));
    }
    
    complexes_.push_back( new RealPairs( pos_real ) );
    
    for(int j=1; j<ix2.size(); ++j) 
    for(int alpha=0; alpha<ix2[j].size(); ++alpha) {
        
        if (ix2[j][alpha]==0) {
            if (length(j)%2) odd_zeros.push_back(length(j));
            else even_zeros.push_back(length(j));
            //result.push_back( new CentralString(length(j)) );
        }
        else if (ix2[j][alpha]>0) {
            int jx2 = guessSum2xJ(big_n, ix2, j, alpha); 
            complexes_.push_back( new String(length(j), jx2) );
        }
    }
    
    if (odd_zeros.size() == 1) {
        // only once central string or origin root
        complexes_.push_back( new CentralString(odd_zeros[0]) );
    }
    else if (odd_zeros.size()==2 && odd_zeros[0]==1 && odd_zeros[1]==3) {
        // 1&3 special kite
        complexes_.push_back( new Kite() );
    }
    else {
        for (int i=0; i<odd_zeros.size(); ++i) {
            // construct generic kite state       
            
            // not implemented
            good_ = false;
            return;
        }
    }
    
    
    if (even_zeros.size() == 1) {
        // only once central string or origin root
        complexes_.push_back( new CentralString(even_zeros[0]) );
    }
    else {
        for (int i=0; i<even_zeros.size(); ++i) {
            // construct generic kite state       
            
            // not implemented
            good_ = false;
            
            return;
        }
    }
    
    return;
}






/* get all roots in complex form 
*/

vector< complex<double> > Solver::getRoots() const
{
 	vector< complex<double> > roots;
    for (int alpha=0; alpha < complexes_.size(); ++alpha) {
        vector< complex<double> > roots_alpha = complexes_[alpha]->getRoots();
        roots.insert(roots.end(), roots_alpha.begin(), roots_alpha.end());
    }
    return roots;
}    


int Solver::numberRoots() const
{
 	int number = 0;
    for (int alpha=0; alpha < complexes_.size(); ++alpha) {
        number += complexes_[alpha]->size();
    }
    return number;
}    


/* calculate the quantum numbers (J) as they are in the complex (non-Takahashi) Bethe equations, from the solution 
*/

vector< complex<double> > Solver::getDirtyJ() const // (const int N, const vector<Roots*>& strings)
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
				qn_J[i] = root_sign*0.25*(big_n_ - (big_n_%4));
			else 
				// even. no zero.
				qn_J[i] = root_sign*0.25*(big_n_ - (2 - big_n_%4));
	    }
	    else {	
		    complex<double> lhs = double(big_n_) * atan(2.*roots[i]);
		    for (int j=0; j<roots.size(); ++j){ 
			    if (j!=i) 
			        lhs -= atan(roots[i] - roots[j]);
		    }
		    qn_J[i] = lhs/PI;
		}	
	}

	return qn_J;
}





