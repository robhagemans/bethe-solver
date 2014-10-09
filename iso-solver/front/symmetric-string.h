#ifndef SYMSTRING_H
#define SYMSTRING_H

#include <vector>
#include <complex>
#include <cmath>

#include "roots.h"
#include "takahashi-string.h"		

class CentralString;


class SymString : public SymRoots {
public:

    SymString(const int big_n, const int string_length, const int sum_jx2, const double lambda = 0.);
    // solving policy, default is used by constructor
    void setPolicy(const double initial_deviation=1e-10, const int steps_no_deviation = 0, const double damping_delta=0.);
    bool setOnAxis(const int n_on_axis);
    
    bool coupleKite(Kite& mate, const int mate_beta);
    
    virtual std::vector<CVal> getRoots() const;
    virtual int size() const;
    
    virtual std::vector<CVal> getComplexQuartets() const;
    virtual std::vector<double> getRealPairs() const; 

    // only nonempty if there are on-axis roots
    virtual std::vector<IVal> getImagPairs() const; 
        
    virtual IterResult iterate(const std::vector<SymRoots*>& all, const int alpha, double& convergence);
    virtual IterResult initiate(const std::vector<SymRoots*>& all, const int alpha);
    virtual void refresh();

    
    friend class CentralString;

private:
    double stepRapidity (const std::vector<SymRoots*>& all, const int alpha);
    bool stepDeviation (const std::vector<SymRoots*>& all, const int alpha, 
                std::vector<double>& new_delta, std::vector<double>& new_epsilon);

    // chain length
    int big_n_;
    
    // number of roots in this string
    int string_length_;          
    
    // sum of bethe quantum numbers for the string, times 2
    int sum_jx2_; 

    bool number_roots_odd_;
    
    // rapidity
    double lambda_;    
    std::vector<double> epsilon_;   
    std::vector<double> delta_;   

    double new_lambda_;    
    std::vector<double> new_epsilon_;   
    std::vector<double> new_delta_;   
    
    // hold a deviation - do not update
    std::vector<bool> hold_;
    
    // number of roots taken to single root on axis
    int n_on_axis_;

    // number of iteration steps passed    
    int iterations_;
    
    Kite* kite_mate_;
    int mate_beta_;
    /* solving policy */

    // size of small deviation to initialize
    double initial_deviation_; 
    // maximum number of steps spent in 'running' phase (exponential growing of lambda to find solvable deviation) 
    //int run_max_; 
    // damping factor for deviance delta
	double damping_delta_; 
    // number of intitial iteration steps in which deviation is skipped  
    int steps_no_deviation_;
};



class CentralString : public SymRoots {
public:
    CentralString(const int big_n, const int string_length);
    // solving policy, default is used by constructor
    void setPolicy(const double initial_deviation=0.001, const int max_simple_step = 0, const double damping_delta=0.);
    // couple with shorter noncentral string - the corresponding roots will be ruled by the mate ann knocked out
    // but play a role in finding roots in this object 
    bool coupleKite(const Kite& mate, const int mate_beta);
    bool coupleString(const SymString& mate, const int mate_beta);
    
    //bool knockOut(const int dummy) {};
    
    // export roots to end user
    virtual std::vector<CVal> getRoots() const;
    // number of roots provided - equals getRoots().size()
    virtual int size() const;
    
    // export roots to other symmetric objects
    virtual std::vector<IVal> getImagPairs() const; 
    //virtual std::vector<double> getImagPairsDelta() const;                         
    //virtual std::vector<int> getImagPairsPos() const; 
    virtual bool hasOrigin() const;
    
    // solving functions
    virtual IterResult iterate(const std::vector<SymRoots*>& all, const int alpha, double& convergence);
    virtual IterResult initiate(const std::vector<SymRoots*>& all, const int alpha);
    virtual void refresh();
    
private:
    bool stepCentralDeviation (const std::vector<SymRoots*>& all, const int alpha, 
                std::vector<double>& new_delta);

    
    void setInitialValues(const double initial_deviation);

    // chain length
    double big_n_;
    
    // number of roots in this string, if nothing knocked out
    int string_length_;    
    
    // sum of bethe quantum numbers for the string, times 2
    int sum_jx2_; 
    
    // number of roots in system is odd?
    // this is set by initiate, once all others are known.
    bool number_roots_odd_;
    
    std::vector<double> delta_;   
    std::vector<double> delta_last_;
    std::vector<double> step_last_;
    // hold the next step before refresh
    std::vector<double> new_delta_;   
       
    
    // hold a deviation - do not update
    std::vector<bool> hold_;
    
    // knockout - do not return roots in the middle of the string
    // number of roots knocked out in the middle of string, at most string_length_-2
    int knockout_;
    const Kite* kite_mate_;        
    const SymString* string_mate_;        
    int mate_beta_;

    // number of iteration steps passed    
    int iterations_;

  
    /* solving policy */

    // size of small deviation to initialize
    double initial_deviation_; 
    // damping factor for deviance delta
	double damping_delta_; 
	// max number of iterations performing simple steps, before using the secant method
	int max_simple_step_;
};


#endif

