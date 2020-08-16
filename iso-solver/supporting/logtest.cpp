#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include <sys/resource.h>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <ieee754.h>

using std::cout;
using std::cerr;
using std::endl;

/** timer **/

class Timer {
private:
    bool system_;
    double start_; 	// initial user time  
    double gettime() const;
public:
	Timer (const bool system_time=false);
	// reset & start
	void start();
	// returns time since stopwatch reset in seconds
	double read() const; 
};


Timer::Timer (const bool system_time): system_(system_time), start_(0)
{	start();  }

void Timer::start() 
{   start_ = gettime();  }

double Timer::read() const
{   return (-start_ + gettime()); }

double Timer::gettime() const
{
    rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    if (system_) 
	    return usage.ru_stime.tv_sec + (usage.ru_stime.tv_usec)/1.0e+6;
	return usage.ru_utime.tv_sec + (usage.ru_utime.tv_usec)/1.0e+6;
}

;

struct LogAdd
{
    double r;
    double inv_r;
    double log_r;
    
    inline LogAdd () : r (1.), log_r(0.), inv_r(1.) {};
    inline void add_log (const double arg)
    {
        r *= arg;
        if (r<1e-10 || r> 1e+10) {
            log_r += log(r);
            r=1.;
        }
    };
    inline void sub_log (const double arg)
    {
        inv_r *= arg;
        if (inv_r<1e-10 || inv_r> 1e+10) {
            log_r -= log(inv_r);
            inv_r=1.;
        }
    };

    inline double read () const { return log_r + log(r) - log(inv_r); };
};


class LogAdd2
{
    static const int limit=IEEE754_DOUBLE_BIAS-100;
    static const int up_limit=IEEE754_DOUBLE_BIAS+100;
    static const double log2_e = M_LOG2E;//1.4426950408889634; // log2(exp(1.)); 

    ieee754_double r;
    ieee754_double inv_r;
    int log2_r;
public:    
    inline LogAdd2 () :  log2_r(0) { r.d=1.; inv_r.d=1.;};
    inline void add_log (const double arg)
    {
        r.d *= arg;
        if (r.ieee.exponent<limit || r.ieee.exponent>up_limit) {
            log2_r += r.ieee.exponent - IEEE754_DOUBLE_BIAS;
            r.ieee.exponent = IEEE754_DOUBLE_BIAS;
        }
    };
    inline void sub_log (const double arg)
    {
        inv_r.d *= arg;
        if (inv_r.ieee.exponent<limit || inv_r.ieee.exponent>up_limit) {
            log2_r -= inv_r.ieee.exponent-IEEE754_DOUBLE_BIAS;
            inv_r.ieee.exponent = IEEE754_DOUBLE_BIAS;
        }
    };

    inline double read () const { return log2_r/log2_e + log(r.d) - log(inv_r.d); };
};


int main()
{
    
    double a;
    double b;
    double c;
    //const int runs = 489910; // here, the last ccyle turns c to inf for the third test, but long double works.
    const int runs =1e+6;    
    Timer timer;
    
    const double rec = 1./RAND_MAX;

    c=0.;        
    srand(0);
    timer.start();
    for (int i=0; i<runs; ++i) {
        a = .5 - rec*rand();
        //b = rec*rand();
        c += log (a*a);
    }
    cout<<c<<" "<<timer.read()<<endl;

    c=0.;        
    srand(0);
    timer.start();
    for (int i=0; i<runs; ++i) {
        a = .5 - rec*rand();
        //b = rec*rand();
        c += 2.*log (abs(a));
    }
    cout<<c<<" "<<timer.read()<<endl;




    c= 0.;
    srand(0);
    timer.start();
    for (int i=0; i<runs; ++i) {
        a = rec*rand();
        b = rec*rand();
        c += log((a*a)/(b*b));
    }
    cout<<c<<" "<<(timer.read())<<endl;
 
    c= 0.;
    srand(0);
    timer.start();
    for (int i=0; i<runs; ++i) {
        a = rec*rand();
        b = rec*rand();
        c += 2.*log(abs(a/b));
    }
    cout<<c<<" "<<(timer.read())<<endl;


exit(0);
/** **/
cout<<endl;
    
    c=0.;        
    srand(0);
    timer.start();
    for (int i=0; i<runs; ++i) {
        a = rec*rand();
        b = rec*rand();
        c += log (a) - log(b);
    }
    cout<<c<<" "<<timer.read()<<endl;
        
    c=0.;        
    srand(0);
    timer.start();
    for (int i=0; i<runs; ++i) {
        a = rec*rand();
        b = rec*rand();
        c += log (a/b);
    }
    cout<<c<<" "<<timer.read()<<endl;
    
    long double d=1.;        
    srand(0);
    timer.start();
    for (int i=0; i<runs; ++i) {
        a = rec*rand();
        b = rec*rand();
        d *= a/b;
    }
    d = log(d);
    cout<<d<<" "<<(timer.read())<<endl;
    

    c=1.;        
    srand(0);
    timer.start();
    for (int i=0; i<runs; ++i) {
        a = rec*rand();
        b = rec*rand();
        c *= a/b;
    }

    d = c;
    d = log(d);
    cout<<d<<" "<<(timer.read())<<endl;
    
    
    d=1.;        
    srand(0);
    timer.start();
    for (int i=0; i<runs; ++i) {
        a = rec*rand();
        b = rec*rand();
        c = a/b;
        d *= c;
    }
    d = log(d);
    cout<<d<<" "<<(timer.read())<<endl;
    
        
        
    c=1.;        
    double e = 0.;
    srand(0);
    timer.start();
    for (int i=0; i<runs; ++i) {
        a = rec*rand();
        b = rec*rand();
        c *= a/b;
        
        if (abs(c)>1e+10 || abs(c)< 1e-10) {
        //if (i%10000) {
            e += log(c);
            c = 1.;
        } 
    }
    e += log(c);
    cout<<e<<" "<<(timer.read())<<endl;
    

    c=1.;        
    e = 0.;
    double df = 1.;
    srand(0);
    timer.start();
    for (int i=0; i<runs; ++i) {
        a = rec*rand();
        b = rec*rand();
        c *= a;
        df *= b;
        if (c< 1e-10 || df< 1e-10) 
        {
            e += log(c) - log(df);
            c = 1.;
            df = 1.;
        }
    }
    e += log(c) - log(df);
    cout<<e<<" "<<(timer.read())<<endl;
    
    
    srand(0);
    LogAdd la;
    LogAdd lb;
    timer.start();
    for (int i=0; i<runs; ++i) {
        a = rec*rand();
        b = rec*rand();
        la.add_log(a);
        lb.add_log(b);
    }
    c= la.read() - lb.read();
    cout<<c<<" "<<(timer.read())<<endl;
    
    
    srand(0);
    LogAdd le;
    timer.start();
    for (int i=0; i<runs; ++i) {
        a = rec*rand();
        b = rec*rand();
        le.add_log(a);
        le.sub_log(b);
    }
    c= le.read();
    cout<<c<<" "<<(timer.read())<<endl;
    
    
    srand(0);
    LogAdd2 lf;
    timer.start();
    for (int i=0; i<runs; ++i) {
        a = rec*rand();
        b = rec*rand();
        lf.add_log(a);
        lf.sub_log(b);
    }
    c= lf.read();
    cout<<c<<" "<<(timer.read())<<endl;
    
/** **/        
    cout<<endl;

    c=0.;        
    srand(0);
    timer.start();
    for (int i=0; i<runs; ++i) {
        a = rec*rand();
        b = rec*rand();
        c += log (a) + log(b);
    }
    cout<<c<<" "<<(timer.read())<<endl;
        
    c=0.;        
    srand(0);
    timer.start();
    for (int i=0; i<runs; ++i) {
        a = rec*rand();
        b = rec*rand();
        c += log (a*b);
    }
    cout<<c<<" "<<(timer.read())<<endl;
        
    d=1.;        
    srand(0);
    timer.start();
    for (int i=0; i<runs; ++i) {
        a = rec*rand();
        b = rec*rand();
        d *= a*b;
    }
    d = log(d);
    
    cout<<d<<" "<<(timer.read())<<endl;
        

    c=1.;        
    srand(0);
    timer.start();
    for (int i=0; i<runs; ++i) {
        a = rec*rand();
        b = rec*rand();
        c *= a*b;
    }
    d = c;
    d = log(d);
    
    cout<<d<<" "<<(timer.read())<<endl;
        


    c=1.;        
    e = 0.;
    srand(0);
    timer.start();
    for (int i=0; i<runs; ++i) {
        a = rec*rand();
        b = rec*rand();
        c *= a*b;
        if (c< 1e-10) {
            e += log(c);
            c = 1.;
        } 
    }
    e += log(c);

    cout<<e<<" "<<(timer.read())<<endl;
       

    LogAdd lc;
    srand(0);
    timer.start();
    for (int i=0; i<runs; ++i) {
        a = rec*rand();
        b = rec*rand();
        lc.add_log(a*b);
    }
    c = lc.read();

    cout<<c<<" "<<(timer.read())<<endl;


    LogAdd2 lz;
    srand(0);
    timer.start();
    for (int i=0; i<runs; ++i) {
        a = rec*rand();
        b = rec*rand();
        lz.add_log(a*b);
    }
    c = lz.read();

    cout<<c<<" "<<(timer.read())<<endl;

/** **/
cout<<endl;
        
    return 0;
}

