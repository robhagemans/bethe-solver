#ifndef SSCAN_H
#define SSCAN_H

#include <vector>
#include "takahashi-string.h"



class SimpleScanner {
public:
    SimpleScanner(const int big_n, const int big_m);

    bool firstState ();
    bool nextState ();

    bool firstBase ();
    bool nextBase ();

    template<class ScanFunc, class SumFunc = NoFunc>
    int scan(ScanFunc per_state, SumFunc per_base = NoFunc(), SumFunc per_scan = NoFunc());

    template<class ScanFunc, class SumFunc>
    int scanBase(ScanFunc per_state, SumFunc per_base, int& not_converged, const std::vector<int>* on_base = 0);


    template<class ScanFunc>
    bool scanState(ScanFunc per_state,  const std::vector<int>* on_base, const int on_state);


private:
    static const int nextNumber (const int number);

    bool firstInType ( const int j, int max_alpha=-1 );

    int big_n_;
    int big_m_;
    std::vector<int> base_;
    std::vector< std::vector<int> > ix2_;
};



template<class ScanFunc>
bool SimpleScanner::scanState(ScanFunc per_state,  const std::vector<int>* on_base, const int on_state)
{
    if (on_base) base_ = *on_base;
    if (!firstState()) return false;

    int num_states = 0;
    do {
        ++num_states;
       	if (num_states==on_state)
       	    return per_state(big_n_, big_m_, base_, ix2_);
    } while (nextState());

    return false;
}


template<class ScanFunc, class SumFunc>
int SimpleScanner::scanBase(ScanFunc per_state, SumFunc per_base, int& not_converged, const std::vector<int>* on_base)
{
    if (on_base) base_ = *on_base;
    if (!firstState()) return 0;

    int not_converged_start = not_converged;
    int num_states = 0;
    do {
        ++num_states;
       	if (!per_state(big_n_, big_m_, base_, ix2_))
       	    ++not_converged;
    } while (nextState());

    per_base(big_n_, big_m_, base_, num_states, not_converged-not_converged_start);
    return num_states;
}


template<class ScanFunc, class SumFunc>
int SimpleScanner::scan(ScanFunc per_state, SumFunc per_base, SumFunc per_scan)
{

    if (!firstBase()) return 0;

    int total_states = 0;
    int not_converged=0;
    do {
        total_states += scanBase(per_state, per_base, not_converged);
    } while (nextBase());

    per_scan(big_n_, big_m_, base_, total_states, not_converged);
    return total_states;
}



#endif
