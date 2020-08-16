
/*
class Kite : public SymRoots {
public:

    Kite(const int big_n, const double level_y=1.); 
    virtual Kite* clone() const;
    
    virtual std::vector< std::complex<double> > getRoots() const;
    
    virtual int size() const;

    virtual std::vector<double> getRealPairs() const;                         
    virtual std::vector<double> getImagPairs() const;                         
    virtual std::vector<double> getImagPairsDelta() const;                         
    virtual std::vector<int> getImagPairsPos() const;                         

    virtual bool iterate(const std::vector<SymRoots*>& all, const int alpha, const double last_convergence, double& convergence);
    virtual bool initiate(const std::vector<SymRoots*>& all, const int alpha);
    virtual void refresh();

    friend class CentralString;
private:
    double big_n_;
    double level_y_;
    
    std::complex<double> step(const std::vector<SymRoots*>& all, const int alpha) const;
    std::complex<double> z_;
    std::complex<double> z_last_;
    std::complex<double> s_last_;
    std::complex<double> z_next_;
};
*/



/** Kite class **/

Kite::Kite (const int big_n, const double level_y) : big_n_(big_n), level_y_(level_y)  { }  

Kite* Kite::clone() const
{   return new Kite(*this);  }

vector< complex<double> > Kite::getRoots() const 
{   return { I*imag(z_), real(z_), -I*imag(z_), -real(z_) };  }

int Kite::size() const
{   return 4;   }

vector<double> Kite::getRealPairs() const
{   return { real(z_) };  }                        

vector<double> Kite::getImagPairs() const                         
{   return { imag(z_) };  }                        

vector<double> Kite::getImagPairsDelta() const                         
{   return { imag(z_) };  }                        

vector<int> Kite::getImagPairsPos() const                         
{   return { 0 };  }                        



complex<double> Kite::step(const vector<SymRoots*>& all, const int alpha) const
{
    double theta = 0.;
    double log_r = 0.;
    int im_log_r_over_pi = 0.;
    
    // scattering with others
    for (int beta=0; beta<all.size(); ++beta) {
        if (beta==alpha) continue;
        
        const bool other_o = all[beta]->hasOrigin();
        const vector<double> other_x = all[beta]->getRealPairs();
        const vector<double> other_y = all[beta]->getImagPairs();
        const vector<complex<double> > other_z = all[beta]->getComplexQuartets();
        
        RealPairs_scatterOthers(real(z_), other_o, other_x, other_y, other_z,  theta);   
        ImagPairs_scatterOthers(imag(z_), other_o, other_x, other_y, other_z,  log_r, im_log_r_over_pi);   
    }    
    
    // now that we have all scattering terms, update our lambdas 
    //return tan( - 0.5*(theta + I*log_r) + 0.5*(big_n_-1.) * ( atan(2.0*real(z_)) + atan(2.0*I*imag(z_)) )  ) ;

cerr<<im_log_r_over_pi<<endl;
    //if (im_log_over_pi%4) return false;

    const double r = sqrt(norm(z_+I))*exp(log_r - (big_n_-1.)*atih(2.0*imag(z_)));
    theta += xi(real(z_), 1.+imag(z_)) -  (big_n_-1.)*atan(2.0*real(z_)) ;     

    return  r*cos(theta) + I*(1.+r*sin(-theta));
}


bool Kite::iterate(const vector<SymRoots*>& all, const int alpha, const double last_convergence, double& convergence)
{
    complex<double> s = step(all, alpha);
    z_next_ = secantStep(z_last_, z_, s_last_, s);
    convergence += norm(z_next_-z_);
    s_last_ = s;
    return true;
}


void Kite::refresh()
{
    z_last_ = z_;
    z_ = z_next_;
}


bool Kite::initiate(const vector<SymRoots*>& all, const int alpha)
{
    z_ = complex<double> (0.01, level_y_+0.001);
    s_last_ = step(all, alpha);    
    z_last_ = z_;
    z_ = complex<double> (0.02, level_y_+0.002);    
    
    return true;
}








