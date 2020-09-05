#ifndef NONSYMSTRING_H
#define NONSYMSTRING_H

#include <vector>
#include <complex>
#include <cmath>

#include "roots.h"
#include "takahashi-string.h"


class NonsymString : public NonsymRoots {
public:
    NonsymString(const int big_n, const int string_length, const int ix2);

    virtual std::vector<CVal> getRoots() const;
    virtual int size() const;

    virtual std::vector<CVal> getComplexPairs() const;
    virtual std::vector<double> getRealRoots() const;

    // rapidities/string centres for each complex
    virtual std::vector<double> getRapidities() const;
    virtual int stringLength() const;


    virtual IterResult iterate(const std::vector<NonsymRoots*>& all, const int alpha, double& convergence);
    virtual IterResult initiate(const std::vector<NonsymRoots*>& all, const int alpha);
    virtual void refresh();

private:
    double stepRapidity (const std::vector<NonsymRoots*>& all, const int alpha);

    bool innerPairIsNarrow( const std::vector<NonsymRoots*>& all, const int alpha) const;

    bool stepDeviation (
                const std::vector<NonsymRoots*>& all, const int alpha,
                std::vector<double>& new_delta, std::vector<double>& new_epsilon);

    void setInitialDeviations(const std::vector<NonsymRoots*>& all, const int alpha);


    // chain length
    int big_n_;

    // number of iteration steps passed
    int iterations_;


    // bethe takahashi quantum number, times 2
    int ix2_;
    // number of roots in this string
    int string_length_;

    // total numbe rof roots in the system, including this string
    int number_roots_;

    // inner sign: is the inner pair of an even string narrow?
    bool inner_narrow_;
    // inner sign has been calculated, use deviations
    bool have_inner_sign_;

    // current
    double lambda_;
    std::vector<double> epsilon_;
    std::vector<double> delta_;

    // next
    double new_lambda_;
    std::vector<double> new_epsilon_;
    std::vector<double> new_delta_;

    // last step value
    double step_lambda_;
    std::vector<double> step_epsilon_;
    std::vector<double> step_delta_;

    // last value
    double last_lambda_;
    std::vector<double> last_epsilon_;
    std::vector<double> last_delta_;


    /* solving policy */
public:
    double damping_dev_;

    // hold a deviation - do not update
    std::vector<bool> hold_;

    // size of small deviation to initialize
    double initial_deviation_;
};




#endif

