#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <cmath>
#include <iostream>
#include "roots.h"
#include "takahashi-string.h"
#include "symmetric-string.h"
#include "nonsymmetric-string.h"


/* interface */





class StateSolver {
public:
    virtual IterResult solve(const int max_iter, const double desired_convergence) = 0;
    
    virtual int iterations() const = 0;
    virtual double convergence() const = 0;

    virtual int size() const = 0;
    virtual std::vector< CVal > getRoots() const = 0;

//    virtual std::vector< std::complex<double> > getRootValues() const = 0;
//    inline std::vector< std::complex<double> > getJs(const int big_n) const { return ::getJs(big_n, getRoots()); };
};


// using function pointer here instead of template argument
// as virtual functions can't be template functions
typedef void (*ReportFunc) (StateSolver*); 
// inline so I can define it here, in the header
inline void noReport(StateSolver* solver) {};

//template<class ReportFunc=NoFunc>
StateSolver* newSolver(
    const int big_n, const std::vector<int>& base, 
    const std::vector< std::vector<int> >& ix2, 
    ReportFunc report=noReport);
    



/* implementation */



// these all refer to takahashi states
// may be worth creating a class for such a thing

bool isSymmetric(const std::vector< std::vector<int> >& ix2);
// faster version where
// vector is assumed to be sorted 0 +1 -1 +2 -2 ...
bool isSymmetricSorted(const std::vector< std::vector<int> >& ix2);
		


template<class RootsType>
class Bunch: public StateSolver {
public:    
    // Bunch takes ownership of the pointers passed and will delete them on destruction (a bit like boost::ptr_vector)
    // - DON'T delete them yourself, best to use new ... in the constructor
    // - DON'T pass &x: if x is statically allocated it will segfault.
    Bunch(std::vector<RootsType*> complexes, ReportFunc report=&noReport): complexes_(complexes), report_(report), iter_(0), conv_(1.) {};
    void setReport(ReportFunc report=&noReport) { report_=report; };
    virtual ~Bunch();
    
    virtual std::vector< CVal > getRoots() const;
        
    virtual int size() const;
    virtual IterResult solve(const int max_iter, const double desired_convergence);

    virtual int iterations() const { return iter_; } ;
    virtual double convergence() const { return conv_; } ;
    
private:                        
    IterResult iterate();
    IterResult initiate();
    void refresh();

private:
    std::vector<RootsType*> complexes_;    
    ReportFunc report_;
    int iter_;
    double conv_;
};



template<class RootsType>
Bunch<RootsType>::~Bunch()
{
    for (int alpha=0; alpha < complexes_.size(); ++alpha) 
        delete complexes_[alpha];   
}

template<class RootsType>
std::vector<CVal> Bunch<RootsType>::getRoots() const
{
 	std::vector<CVal> roots;
    for (int alpha=0; alpha < complexes_.size(); ++alpha) {
        std::vector<CVal> roots_alpha = complexes_[alpha]->getRoots();
        roots.insert(roots.end(), roots_alpha.begin(), roots_alpha.end());
    }
    return roots;
}    


template<class RootsType>
int Bunch<RootsType>::size() const
{
 	int number = 0;
    for (int alpha=0; alpha < complexes_.size(); ++alpha) {
        number += complexes_[alpha]->size();
    }
    return number;
}    

template<class RootsType>
IterResult Bunch<RootsType>::initiate() 
{
    IterResult result = IterResult::iter_ok;
    for (int alpha=0; alpha < complexes_.size(); ++alpha) {
        IterResult error = complexes_[alpha]->initiate(complexes_, alpha);
        // highest error value dominates
        if (error>result) result = error; 
    }
    return result;
}

template<class RootsType>
IterResult Bunch<RootsType>::iterate()
{
    IterResult result = IterResult::iter_ok;
    for (int alpha=0; alpha < complexes_.size(); ++alpha) {
        IterResult error = complexes_[alpha]->iterate(complexes_, alpha, conv_); 
        // highest error value dominates
        if (error>result) result = error; 

    }
    return result;
}


template<class RootsType>
void Bunch<RootsType>::refresh()
{
    for (int alpha=0; alpha < complexes_.size(); ++alpha) {
        complexes_[alpha]->refresh(); 
    }
}

template<class RootsType>
IterResult Bunch<RootsType>::solve(const int max_iter, const double desired_convergence) 
{
    IterResult iter_result = IterResult::iter_err;
    
    if (0==size()) return iter_result;
    if (fatal(iter_result = initiate())) 
        return iter_result;

    report_(this);
    
    for (; iter_<max_iter; ++iter_) {
        
        conv_ = 0.;
        // calculate next value
        if (fatal(iter_result = iterate())) 
            return iter_result;
        
        // replace old value with new
        refresh();
        
        report_(this);
        
        if (!finite(conv_)) 
            return IterResult::iter_err_nonfinite;        
            
        if (conv_ < desired_convergence) 
            return iter_result;        
    }
    
    // not converged if here
    return IterResult::solver_not_converged;
}



#endif

