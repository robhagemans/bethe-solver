#ifndef SYMSTRING_H
#define SYMSTRING_H

#include <vector>
#include <complex>
#include <cmath>

#include "roots.h"
#include "takahashi-string.h"		

class CentralString;

class String : public SymRoots {
public:

    String(const int big_n, const int string_length, const int sum_jx2, const double lambda = 0.);
    // solving policy, default is used by constructor
    void setPolicy(const double initial_deviation=1e-10, const int steps_no_deviation = 0, const double damping_delta=0.);
    
    virtual std::vector< std::complex<double> > getRoots() const;
    virtual int size() const;
    
    virtual std::vector< std::complex<double> > getComplexQuartets() const;
    virtual std::vector<double> getRealPairs() const; 

        
    virtual bool iterate(const std::vector<SymRoots*>& all, const int alpha, const double last_convergence, double& convergence);
    virtual bool initiate(const std::vector<SymRoots*>& all, const int alpha);
    virtual void refresh();

    
    static long double thetaOthers(const double x, const double y, 
                            bool has_origin, const std::vector<double>& other_x, 
                            const std::vector<double>& other_y, const std::vector<std::complex<double> >& other_z);
                            
    static long double logRSqOthers(const double x, const double y, 
                            bool has_origin, const std::vector<double>& other_x, 
                            const std::vector<double>& other_y, const std::vector<std::complex<double> >& other_z);
                            

    friend class CentralString;

private:
    double stepRapidity (const std::vector<SymRoots*>& all, const int alpha);
    bool stepDeviation (const std::vector<SymRoots*>& all, const int alpha, 
                std::vector<double>& new_delta, std::vector<double>& new_epsilon);

    
    // chain length
    double big_n_;
    
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


    // number of iteration steps passed    
    int iterations_;
    
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
    bool coupleKite(const Kite& mate, const int mate_beta_);
    bool coupleString(const String& mate, const int mate_beta_);
    
    bool knockOut(const int dummy) {};
    
    // export roots to end user
    virtual std::vector< std::complex<double> > getRoots() const;
    // number of roots provided - equals getRoots().size()
    virtual int size() const;
    
    // export roots to other symmetric objects
    virtual std::vector<double> getImagPairs() const; 
    virtual std::vector<double> getImagPairsDelta() const;                         
    virtual std::vector<int> getImagPairsPos() const; 
    virtual bool hasOrigin() const;
    
    // solving functions
    virtual bool iterate(const std::vector<SymRoots*>& all, const int alpha, const double last_convergence, double& convergence);
    virtual bool initiate(const std::vector<SymRoots*>& all, const int alpha);
    virtual void refresh();
    
    // get Bethe quantum numbers
    virtual std::vector<std::complex<double> > getCleanerJ (const std::vector<SymRoots*>& all, const int alpha) const;
    virtual std::vector<int> getCleanJx2 (const std::vector<SymRoots*>& all, const int alpha, double& sq_error) const;

private:
    bool stepCentralDeviation (const std::vector<SymRoots*>& all, const int alpha, 
                std::vector<double>& new_delta);

    void scatterOthers(const std::vector<SymRoots*>& all, int& theta, long double& log_rsq, const int alpha, const int a) const;
    
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
    
    // deviance
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
    const String* string_mate_;        
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



class KiteString : public SymRoots {
public:

    KiteString(const int big_n, const int string_length);// const int sum_jx2 );
    // solving policy, default is used by constructor
    void setPolicy(const double initial_deviation=1e-10, const double damping_delta=0.);
    
    virtual std::vector< std::complex<double> > getRoots() const;
    virtual int size() const;
    
    virtual std::vector< std::complex<double> > getComplexQuartets() const;
    virtual std::vector<double> getRealPairs() const; 
    virtual std::vector<double> getImagPairs() const; 

        
    virtual bool iterate(const std::vector<SymRoots*>& all, const int alpha, const double last_convergence, double& convergence);
    virtual bool initiate(const std::vector<SymRoots*>& all, const int alpha);
    virtual void refresh();


private:
    double stepRapidity (const std::vector<SymRoots*>& all, const int alpha);

    bool stepDeviation (const std::vector<SymRoots*>& all, const int alpha, 
                        std::vector<double>& new_delta, std::vector<double>& new_epsilon);

    void setInitialValues(const double initial_eps, const double initial_del);

    // chain length
    double big_n_;
    
    // number of roots in this string
    int string_length_;          
    
    // sum of bethe quantum numbers for the string, times 2
    //int sum_jx2_; 

    bool number_roots_odd_;
    
    std::vector<double> epsilon_;   
    std::vector<double> delta_;   

    std::vector<double> new_epsilon_;   
    std::vector<double> new_delta_;   
    
    std::vector<double> last_epsilon_;   
    std::vector<double> last_delta_;   
    
    std::vector<double> step_epsilon_;   
    std::vector<double> step_delta_;   
    
    
    // hold a deviation - do not update
    std::vector<bool> hold_;


    // number of iteration steps passed    
    int iterations_;
    
    /* solving policy */

    // size of small deviation to initialize
    double initial_deviation_;
    // maximum number of steps spent in 'running' phase (exponential growing of lambda to find solvable deviation) 
    //int run_max_; 
    // damping factor for deviance delta
	double damping_delta_; 
};


class String2 : public SymRoots {
public:

    String2(const int big_n, const int string_length, const int sum_jx2, const double lambda = 0.);
    // solving policy, default is used by constructor
    void setPolicy(const double initial_deviation=1e-10, const int steps_no_deviation = 0, const double damping_delta=0.);
    bool setOnAxis(const int n_on_axis);
    
    virtual std::vector< std::complex<double> > getRoots() const;
    virtual int size() const;
    
    virtual std::vector< std::complex<double> > getComplexQuartets() const;
    virtual std::vector<double> getRealPairs() const; 
    virtual std::vector<double> getImagPairs() const; 

        
    virtual bool iterate(const std::vector<SymRoots*>& all, const int alpha, const double last_convergence, double& convergence);
    virtual bool initiate(const std::vector<SymRoots*>& all, const int alpha);
    virtual void refresh();

    
    static long double thetaOthers(const double x, const double y, 
                            bool has_origin, const std::vector<double>& other_x, 
                            const std::vector<double>& other_y, const std::vector<std::complex<double> >& other_z);
                            
    static long double logRSqOthers(const double x, const double y, 
                            bool has_origin, const std::vector<double>& other_x, 
                            const std::vector<double>& other_y, const std::vector<std::complex<double> >& other_z);
                            
private:
    double stepRapidity (const std::vector<SymRoots*>& all, const int alpha);
    bool stepDeviation (const std::vector<SymRoots*>& all, const int alpha, 
                std::vector<double>& new_delta, std::vector<double>& new_epsilon);

    
    // chain length
    double big_n_;
    
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


#endif

