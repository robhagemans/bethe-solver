#ifndef ROOTS_H
#define ROOTS_H

#include <complex>
#include <vector>
#include <cmath>

#include "basics.h"

// would like to  un-implement the abs in global namespace (as non-std abs turns double to ints, hard to debug)
// at least this way things work as expected
using std::abs;



/** root value structs **/

struct IVal 
{
    int pos;
    double del;
    
    inline double imag() const { return 0.5*pos+del; }; 
};

struct CVal 
{
    double lam;
    double eps;
    int pos;
    double del;
    
    IVal ival() const { return {pos, del }; };
    
    inline double real() const { return lam+eps; }; 
    inline double imag() const { return 0.5*pos+del; };  
    inline std::complex<double> val() const { return std::complex<double>(lam+eps, 0.5*pos+del); }; 
};


// we use c++ convention that norm is the square of the absolute value
inline double norm(const CVal& a)
{   return sq(a.real()) + sq(a.imag()); }

inline CVal sum(const CVal& a, const CVal& b)
{   return { a.lam+b.lam, a.eps+b.eps, a.pos+b.pos, a.del+b.del }; }

inline IVal sum(const IVal& a, const IVal& b)
{   return { a.pos+b.pos, a.del+b.del }; }

inline CVal sum(const CVal& a, const double b)
{   return { a.lam+b, a.eps, a.pos, a.del }; }
inline CVal sum(const double a, const CVal& b)
{   return { a+b.lam, b.eps, b.pos, b.del }; }

inline CVal sum(const CVal& a, const IVal& b)
{   return { a.lam, a.eps, a.pos+b.pos, a.del+b.del }; }
inline CVal sum(const IVal& a, const CVal& b)
{   return sum (b, a);  }

inline CVal sum(const double a, const IVal& b)
{   return { a, 0., b.pos, b.del }; }
inline CVal sum(const IVal& a, const double b)
{   return sum(b, a); }



inline CVal dif(const CVal& a, const CVal& b)
{   return { a.lam-b.lam, a.eps-b.eps, a.pos-b.pos, a.del-b.del }; }

inline IVal dif(const IVal& a, const IVal& b)
{   return { a.pos-b.pos, a.del-b.del }; }

inline CVal dif(const CVal& a, const IVal& b)
{   return { a.lam, a.eps, a.pos-b.pos, a.del-b.del }; }
inline CVal dif(const IVal& a, const CVal& b)
{   return { -b.lam, -b.eps, a.pos-b.pos, a.del-b.del }; }

inline CVal dif(const CVal& a, const double b)
{   return { a.lam-b, a.eps, a.pos, a.del }; }
inline CVal dif(const double a, const CVal& b)
{   return { a-b.lam, -b.eps, -b.pos, -b.del }; }

inline CVal dif(const double a, const IVal& b)
{   return { a, 0., -b.pos, -b.del }; }
inline CVal dif(const IVal& a, const double b)
{   return { -b, 0., a.pos, a.del }; }



inline CVal mul(const CVal& a, const int n)
{   return { a.lam*n, a.eps*n, a.pos*n, a.del*n }; }
inline IVal mul(const IVal& a, const int n)
{   return { a.pos*n, a.del*n }; }


// returns y + I*n
inline CVal add_n_i(const CVal& z, int n)
{   return { z.lam, z.eps, z.pos+2*n, z.del };  } 
inline IVal add_n_i(const IVal& y, int n)
{   return { y.pos+2*n, y.del };  } 


inline CVal neg(const CVal& a)
{   return { -a.lam, -a.eps, -a.pos, -a.del }; }
inline CVal conj(const CVal& a)
{   return { a.lam, a.eps, -a.pos, -a.del }; }

inline IVal neg(const IVal& y)
{   return { -y.pos, -y.del };  } 
inline IVal conj(const IVal& y)
{   return neg(y);  } 


/** HPIVal **/
// number of half-pi factors
// for use in arguments

struct HPIVal 
{
    int n; 
    double eps;
    
    inline double real() const { return M_PI_2*(n+eps); }; 
    inline double val() const { return real(); };
    
    inline HPIVal& add(const HPIVal& b) { n+=b.n; eps+=b.eps; return *this; };
    inline HPIVal& sub(const HPIVal& b) { n-=b.n; eps-=b.eps; return *this; };

    inline HPIVal& add_half_pi(const int b) { n+=b; return *this; };
    inline HPIVal& sub_half_pi(const int b) { n-=b; return *this; };

};

inline HPIVal sum(const HPIVal& a, const HPIVal& b)
{   return { a.n+b.n, a.eps+b.eps }; }

inline HPIVal mul(const HPIVal& a, const int n)
{   return { a.n*n, a.eps*n }; }

inline HPIVal div(const HPIVal& a, const int n)
{   return { a.n/n, (a.eps +(a.n%n))/double(n) }; }

inline HPIVal neg(const HPIVal& a)
{   return { -a.n, -a.eps }; }

inline double sin(const HPIVal& a) 
{ 
     return a.n%2==0 
        ? ( a.n%4==0 ? std::sin(a.eps*M_PI_2) : -std::sin(a.eps*M_PI_2) ) 
        : ( a.n%4==1 ? std::cos(a.eps*M_PI_2) : -std::cos(a.eps*M_PI_2) );
}

inline double cos(const HPIVal& a) 
{ 
     return a.n%2==0 
        ? ( a.n%4==0 ? std::cos(a.eps*M_PI_2) : -std::cos(a.eps*M_PI_2) ) 
        : ( a.n%4==1 ? -std::sin(a.eps*M_PI_2) : std::sin(a.eps*M_PI_2) );
}

inline double tan(const HPIVal& a) 
{ 
     return (a.n%2) ? (-1./std::tan(a.eps*M_PI_2)) : (std::tan(a.eps*M_PI_2));
}


/** roots and Bethe equations **/

// roots
std::vector< std::complex<double> > getValues(const std::vector<CVal>& roots);

// bethe quantum numbers
std::vector< std::complex<double> > getJs(const int big_n, const std::vector<CVal>& roots);

std::vector< std::complex<double> > getJs_BT(const int big_n, const std::vector<CVal>& roots);


// error is given as average of square error per quantum number
std::vector<int> getJx2 (const int big_n, const std::vector<CVal>& roots, double* error_per_number=0, bool* parity_error=0) ;

// bethe quantum numbers, sanity check
// some cancellations will not work out accurate so you'll end up with imaginary parts if deltas are small
std::vector< std::complex<double> > getDirtyJs(const int big_n, const std::vector<CVal>& roots) ;




/* modular arithmetic */

// double number moduulo 
template<typename number>
inline number modulo (number value, const number period) 
{ 
	while (value < 0) value += period; 
	while (value >= period) value -= period; 
	return value; 
}

// mod arithmetic: not quite the same as % for negative values
// int specialisation of template 
template<>
inline int modulo<int> (int value, const int period) 
{ 	return (value %= period)>=0 ? value : value+period; }
 



inline bool zero(const double x)
{  return x==0.; }
//{  return std::abs(x)<EPSILON_ATAN; }

inline bool zero_re(const CVal& z)
{    return zero(z.lam) && zero(z.eps); }

inline bool zero_re(const CVal& z, const CVal& z_k)
{    return zero_re(dif(z, z_k)); }

inline bool zero_im(const CVal& z)
{    return z.pos==0 && zero(z.del);  }

inline bool zero_im(const CVal& z, const CVal& z_k)
{    return zero_im(dif(z, z_k));  }

inline bool zero_im(const IVal& z)
{    return z.pos==0 && zero(z.del);  }


inline bool abs_greater_than_one (const IVal& y)
{    return ((std::abs(y.pos)==2) ? ((y.pos>0)?(y.del>0.):(y.del<0.)) : (std::abs(y.imag())>1.)); }

inline bool abs_equals_one (const IVal& y)
{    return ((std::abs(y.pos)==2) && (y.del==0.)); }



/** scattering **/


// THIS IS NOT IMPLEMENTED AND WILL GIVE A LINKER ERROR IF USED 
// this is intentional: it avoids inadvertent use of std::atan through a using directive
// std::atan can still be used if namespace explicitly specified
std::complex<double> atan(const std::complex<double> argument);




/** scattering: real parts **/

// this is just the atan as the imag part is zero
// however, this version separates the half-pi term
inline HPIVal re_atan (const double x)
{   return { (std::abs(x)<1 ? 0 : isgn(x)), M_1_PI*atan( 2*x /(1-sq(x)) )  };  }


// should be NAN if abs_equals_one(y)
// but ints can't actually hold NAN
inline int re_atan_x2_div_pi(const IVal& y)
{   return (abs_greater_than_one(y) ? isgn(y.imag()) : 0); }

inline HPIVal re_atan(const IVal& y)
{   return abs_equals_one(y) ? HPIVal { 0, NAN } : HPIVal { re_atan_x2_div_pi(y), 0. }; }


//inline double re_atan_x2(const CVal& z)
//{   return xi( z.real(), 0.5*(2+z.pos)+z.del ) +  xi( z.real(), 0.5*(2-z.pos)-z.del )
//            + ( (zero_re(z) && !zero_im(z)) ? re_atan_x2(z.ival()) : 0.);   }




inline HPIVal re_atan (const CVal& z)
{   return zero_re(z) 
                ? ( re_atan(z.ival()) )
                : ( HPIVal { 
                        ( norm(z)<1 ? 0 : isgn(z.real()) ),  
                        M_1_PI*atan( 2.*z.real() / (-add_n_i(z.ival(),-1).imag() * add_n_i(z.ival(),1).imag() - sq(z.real())) ) 
                    } )
                ;
}
    


/** real half terms (xi) **/

// atan(a+ib) + atan(a-ib) == xi(a, 1+b) + xi(a,1-b)
// re(atan a+ib) = atan(a+ib) + atan(a-ib)   except on the imag axis
//
// old definition:
//   double xi (const double epsilon, const double delta)  { return (atan(epsilon/delta) + ((delta<0.)?M_PI*sgn(epsilon):0)); }	

// argument, defined as im log(i*z)
// branch cut along negative real axis, with the negative re axis itself having an argumnet of +pi,
// limit from above +pi, and the limit from below is -pi.
// returned as an HPIVal such that the eps value is always less than 0.5pi absolute, and the n value holds all terms 0.5pi.   
inline HPIVal im_log_i (const CVal& z)
{
    return zero_re(z) && zero_im(z)  
        ? ( HPIVal { 0, NAN } )  
        : ( std::abs(z.real()) < std::abs(z.imag()) )
            ? (  
                HPIVal { (z.imag()<0 ? 0: 2*isgn_plus(z.real()) ) , -M_2_PI*atan(z.real()/z.imag()) } 
              )
            : (
                HPIVal { isgn(z.real()), M_2_PI*atan(z.imag()/z.real()) } 
              );        
}




// xi+(z) = xi(x, 1+y) = - im log (1-iz) = -im log i(-i-z) except on the imag axis where xi+(iy) = 0
inline HPIVal xi_plus (const CVal& z)
{
    return zero_re(z)
        ? ( HPIVal { 0, (z.pos==-2 && z.del==0 ? NAN : 0) } ) 
        : ( neg(im_log_i(neg(add_n_i(z, 1)))) );   
}

// xi-(z) = xi(x, 1-y) = + im log (1+iz) =im log i (-i+z) except on the imag axis where xi-(iy) = 0
inline HPIVal xi_minus (const CVal& z)
{
    return zero_re(z) 
        ? ( HPIVal { 0, (z.pos==2 && z.del==0 ? NAN : 0) } ) 
        : ( im_log_i(add_n_i(z, -1)) );   
}



/** scattering: imaginary parts **/


//inline double im_atan_x4(const IVal& y)
//{   return (zero_im(y))? 0. : log ( ( sq(0.5*(2+y.pos)+y.del) ) / ( sq(0.5*(2-y.pos)-y.del) )   );  }

inline double im_atan_x4(const IVal& y)
{   
    return zero_im(y) 
        ? ( 0. ) 
        : ( y.pos==0 && std::abs(y.del)<1. 
         ? ( 4.*atanh(y.del) )     // better precision if y.pos==0
         : ( 2.*log(std::abs( (0.5*(2+y.pos)+y.del) / (0.5*(2-y.pos)-y.del) )) )
        );  
}

// what if re(z) is close to one?
inline double im_atan_x4(const CVal& z)
{   
    return zero_re(z) 
        ? ( im_atan_x4(z.ival()) )
        : ( zero_im(z) ? 0 :
                log ( ( sq(z.real()) + sq(0.5*(2+z.pos)+z.del) ) / 
                      ( sq(z.real()) + sq(0.5*(2-z.pos)-z.del) )   )   
          );  
            
}     
     
     
/** imag half terms (zeta) **/

// atan(a+ib) - atan(a-ib) == i (zeta(a,1+b) - zeta(a,1-b))
//
// old definition:
//  inline double zeta(const double epsilon, const double delta)   {   return 0.5*log( epsilon*epsilon + delta*delta );  }


// 2*zeta(0, y+1)
inline double zeta_plus_x2(const IVal& y)
{   
    return (y.pos==0 && std::abs(y.del)<1.)
            ? ( 
                2.*log1p(y.del)
              )  
            : (
                2.*log(std::abs( 0.5*(2+y.pos)+y.del ))
              );  
}  

// 2*zeta(0, y-1)
inline double zeta_minus_x2(const IVal& y)
{   
    return (y.pos==0 && std::abs(y.del)<1.)
            ? ( 
                2.*log1p(-y.del)
              )
            : (  
                2.*log(std::abs( 0.5*(2-y.pos)-y.del ))
              );  
}  

// 2*zeta(x, y+1)
inline double zeta_plus_x2(const CVal& z)
{   
    return zero_re(z) 
            ? ( zeta_plus_x2(z.ival()) )
            : ( log(sq(z.real()) + sq(0.5*(2+z.pos)+z.del)) );  
}  

// 2*zeta(x, y-1)
inline double zeta_minus_x2(const CVal& z)
{   return zero_re(z) 
            ? ( zeta_minus_x2(z.ival()) )
            : ( log(sq(z.real()) + sq(0.5*(2-z.pos)-z.del)) );  
}  




/** scattering sums **/

// symmetric

HPIVal sum_re_atan(const double x,  
        const bool has_origin, const std::vector<double>& other_x, 
        const std::vector<IVal>& other_y, const std::vector<CVal>& other_z);

double sum_im_atan(const IVal y,  
        const bool has_origin, const std::vector<double>& other_x, 
        const std::vector<IVal>& other_y, const std::vector<CVal>& other_z, 
        int& adds_re_atan_x4_div_pi); 

HPIVal sum_re_atan(const CVal& z, 
        const bool has_origin, const std::vector<double>& other_x, 
        const std::vector<IVal>& other_y, const std::vector<CVal>& other_z) ;

double sum_im_atan_x4(const CVal& z, 
        const bool has_origin, const std::vector<double>& other_x, 
        const std::vector<IVal>& other_y, const std::vector<CVal>& other_z) ;
                      
// nonsymmetric

HPIVal sum_re_atan(const double x, const std::vector<double>& other_x, const std::vector<CVal>& other_z) ;
HPIVal sum_re_atan(const CVal& z, const std::vector<double>& other_x, const std::vector<CVal>& other_z) ;
double sum_im_atan_x4(const CVal& z, const std::vector<double>& other_x, const std::vector<CVal>& other_z) ;





/** iteration helpers **/

template <typename number>
number secantStep(
        const number z_last, const number z, 
        const number s_last, const number s)
{ 
        return ((number)0 == z_last -s_last -z +s) 
                        ? s   
                        : ( z*(z_last-s_last) - z_last*(z-s) ) /   ( z_last -s_last -z +s );
}

template <typename number>
number dampedStep(
        const number old, 
        const number step, 
        const double damping)
{ 
        return (1.-damping)*step + damping*old;
}



/** solving result **/

// iteration result
enum class IterResult {
    iter_ok,                // iteration went normal
// nonfatal:    
    iter_constrained,       // value was constrained, not a valid solution but ok to keep iterating
    iter_takahashi,         // value was obtained by bethe-takahashi approximation, ok to keep iterating
    iter_err_sign,          // sign error - LHS = -RHS, sometimes ok to keep iterating
// fatal:    
    iter_err_neg_norm,      // negative argument to square root, occasionally ok to keep iterating if implemented with abs
    iter_err_nonfinite,     // non-finite result
    iter_err_diverged,      // divergence cutoff
    iter_err,               // unspecified error
    solver_not_converged    // solver level only - no convergence but no other errors
};

// result is fatal error for Solver
inline bool fatal(IterResult iter_result)
{   
    return !(  iter_result==IterResult::iter_ok
                || iter_result==IterResult::iter_err_sign 
                || iter_result==IterResult::iter_constrained 
                || iter_result==IterResult::iter_takahashi );
}  


inline const char* message(IterResult iter_result)
{
    switch(iter_result) {
        case IterResult::iter_ok: return "ok";
        case IterResult::iter_constrained: return "active constraint";
        case IterResult::iter_takahashi: return "Bethe-Takahashi approximation";
        case IterResult::iter_err_nonfinite: return "non-finite value";
        case IterResult::iter_err_diverged: return "diverged";
        case IterResult::iter_err_neg_norm: return "negative argument to sqrt norm";
        case IterResult::iter_err_sign: return "sign error";
        case IterResult::iter_err: return "unspecified error";
        case IterResult::solver_not_converged: return "not converged";
    }
    return "unknown error value";
}



/** Roots ABC **/

class Roots {
public:
    // these would be the same for non-symmetric solvers, so single-count every root    
    virtual std::vector<CVal> getRoots() const = 0;
    virtual int size() const = 0;

    // commit changes of last iteration.
    virtual void refresh() { };
};





class NonsymRoots: public Roots {
public:
    inline virtual std::vector<double> getRealRoots() const { return {}; };                         
    inline virtual std::vector<CVal> getComplexPairs() const  { return {}; };                         
    
    virtual IterResult iterate(const std::vector<NonsymRoots*>& all, const int alpha_self, double& convergence) = 0;
    virtual IterResult initiate(const std::vector<NonsymRoots*>& all, const int alpha_self) = 0;

    // rapidities/string centres for each complex 
    virtual std::vector<double> getRapidities() const = 0; 
    virtual int stringLength() const = 0;    
};


class SymRoots: public Roots {
public:
    inline virtual std::vector<double> getRealPairs() const { return {}; };                         
    inline virtual std::vector<IVal> getImagPairs() const { return {}; };                         
    inline virtual std::vector<CVal> getComplexQuartets() const { return {}; };                         
    inline virtual bool hasOrigin() const { return false; };                         
    
    virtual IterResult iterate(const std::vector<SymRoots*>& all, const int j_self,  double& convergence) = 0;
    virtual IterResult initiate(const std::vector<SymRoots*>& all, const int j_self) = 0;
};



/** special cases **/

class OriginRoot: public SymRoots {
public:    
    virtual OriginRoot* clone() const {  return new OriginRoot; };
    virtual std::vector<CVal> getRoots() const { return { {0.,0.,0,0.} }; };
    virtual int size() const { return 1; };
    virtual bool hasOrigin() const { return true; };                         
    virtual IterResult iterate(const std::vector<SymRoots*>& all, const int j_self, double& convergence) 
            { return IterResult::iter_ok; };
    virtual IterResult initiate(const std::vector<SymRoots*>& all, const int j_self) { return IterResult::iter_ok; };
};

        
class OriginPair: public SymRoots {
    virtual OriginPair* clone() const { return new OriginPair; };
    virtual std::vector<CVal> getRoots() const { return { {0.,0.,1,0.}, {0.,0.,-1,0.} };  };
    virtual int size() const { return 2; };
    virtual std::vector<IVal> getImagPairs() const { return { {1, 0.} }; };
    
    virtual IterResult iterate(const std::vector<SymRoots*>& all, const int j_self, double& convergence)
             { return IterResult::iter_ok; };
    virtual IterResult initiate(const std::vector<SymRoots*>& all, const int j_self) { return IterResult::iter_ok; };
};



class RealPairs : public SymRoots {
public:

    RealPairs(const int big_n, const std::vector<int>& jx2, const std::vector<double>& lambda);
    RealPairs(const int big_n, const std::vector<int>& jx2);
    
    virtual RealPairs* clone() const;

    virtual std::vector<CVal> getRoots() const;
    virtual int size() const;
 
    virtual std::vector<double> getRealPairs() const;                         
    
    virtual IterResult iterate(const std::vector<SymRoots*>& all, const int alpha, double& convergence);
    virtual IterResult initiate(const std::vector<SymRoots*>& all, const int alpha);
    virtual void refresh();
    
private:
    double step(const std::vector<SymRoots*>& all, const int alpha, const int i) const;

    double big_n_;
    std::vector<int> jx2_;           // 2 times the bethe quantum number
    std::vector<double> lambda_;     // current value
    std::vector<double> lambda_last_; // previous value 
    std::vector<double> s_last_;    // previous value of step function
    std::vector<double> z_next_;    // next lambda
};


class NonsymRealRoots : public NonsymRoots {
public:

    NonsymRealRoots(const int big_n, const std::vector<int>& jx2, const std::vector<double>& lambda);
    NonsymRealRoots(const int big_n, const std::vector<int>& jx2);

    // add a root to the vector
    virtual void add(const int jx2, const double lambda);
    
    virtual std::vector<CVal> getRoots() const;
    virtual int size() const;
 
    virtual std::vector<double> getRealRoots() const;                         

    // takahashi string hypothesis values
    virtual std::vector<double> getRapidities() const; 
    virtual int stringLength() const;
  
    virtual IterResult iterate(const std::vector<NonsymRoots*>& all, const int alpha, double& convergence);
    virtual IterResult initiate(const std::vector<NonsymRoots*>& all, const int alpha);
    virtual void refresh();
    
private:
    double big_n_;
    double step(const std::vector<NonsymRoots*>& all, const int alpha, const int i) const;

    std::vector<int> jx2_;              // 2 times the bethe quantum number
    std::vector<double> lambda_;        // current value
    std::vector<double> lambda_last_;   // previous value 
    std::vector<double> s_last_;        // previous value of step function
    std::vector<double> z_next_;        
};




class ImagPairs : public SymRoots {
public:

    ImagPairs(const int big_n, const int dummy, const std::vector<int>& jx2);
    ImagPairs(const int big_n, const std::vector<double>& level);
    
    virtual ImagPairs* clone() const;

    virtual std::vector<CVal> getRoots() const;
    virtual int size() const;
    
    virtual std::vector<IVal> getImagPairs() const;                         
                         
    virtual IterResult iterate(const std::vector<SymRoots*>& all, const int alpha, double& convergence);
    virtual IterResult initiate(const std::vector<SymRoots*>& all, const int alpha);
    virtual void refresh();

//  currently refreshes immediately on step - ie before others do their step
//  virtual void refresh();

private:

    double step(const std::vector<SymRoots*>& all, const int alpha, const int i) const;

    double big_n_;
    std::vector<double> level_;   // string level/ minimum imaginary level for this root
    
    std::vector<double> im_lambda_;
    std::vector<double> im_lambda_last_; // previous value 
    std::vector<double> s_last_;   // previous value of step function
    std::vector<double> z_next_;   // previous value of step function
};



class ImagPairs2 : public SymRoots {
public:

    ImagPairs2(const int big_n, const std::vector<int>& pos);
    
    virtual std::vector<CVal> getRoots() const;
    virtual int size() const;
    
    virtual std::vector<IVal> getImagPairs() const;                         
                         
    virtual IterResult iterate(const std::vector<SymRoots*>& all, const int alpha, double& convergence);
    virtual IterResult initiate(const std::vector<SymRoots*>& all, const int alpha);
    virtual void refresh();
        
private:

    double step(const std::vector<SymRoots*>& all, const int alpha, const int i) const;

    double big_n_;
    std::vector<int> pos_;   
    
    std::vector<double> delta_;
    std::vector<double> delta_next_; // previous value 
    std::vector<double> delta_last_; // previous value 
    std::vector<double> s_last_;   // previous value of step function
};



class CentralString;
class SymString;


class Kite : public SymRoots {
public:

    Kite (const int big_n, const int pos);
    void setInitial(const double x, const double del, const double x1=0., const double del1=0.);

    virtual std::vector<CVal> getRoots() const;
    
    virtual int size() const;

    virtual std::vector<double> getRealPairs() const;                         
    virtual std::vector<IVal> getImagPairs() const;                           

    virtual IterResult iterate(const std::vector<SymRoots*>& all, const int alpha, double& convergence);
    virtual IterResult initiate(const std::vector<SymRoots*>& all, const int alpha);
    virtual void refresh();

    
    friend class CentralString;
    friend class SymString;
private:
    
    bool step(const std::vector<SymRoots*>& all, const int alpha, double& new_x, double& new_y) const;

    double big_n_;
    
    double x_;
    int pos_;
    double del_;

    std::complex<double> z_last_;
    std::complex<double> s_last_;
    std::complex<double> z_next_;
};





#endif

