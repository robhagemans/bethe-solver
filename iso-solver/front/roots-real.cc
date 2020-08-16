
#include "roots.h"

using namespace std;


/* RealPairs class */

RealPairs::RealPairs(const int big_n, const std::vector<int>& jx2)
      : big_n_(big_n), jx2_(jx2), lambda_(jx2.size()), lambda_last_(jx2.size()), s_last_(jx2.size()), z_next_(jx2.size())
{ }

RealPairs::RealPairs(const int big_n, const vector<int>& jx2, const vector<double>& lambda)
    :   big_n_(big_n), jx2_(jx2), lambda_(lambda), lambda_last_(jx2.size()), s_last_(jx2.size()), z_next_(jx2.size())
{ }

RealPairs* RealPairs::clone() const
{   return new RealPairs(*this); }

int RealPairs::size() const
{   return lambda_.size()*2; }

vector<double> RealPairs::getRealPairs() const
{   return lambda_;  }


vector<CVal> RealPairs::getRoots() const
{
    vector<CVal> roots (lambda_.size()*2);
    for(int i = 0; i<lambda_.size(); ++i) {
        roots[i] = { lambda_[i], 0.,0,0. };
        roots[2*lambda_.size()-i-1] = { -lambda_[i], 0.,0,0. };
    }
    return roots;
}



double RealPairs::step(const vector<SymRoots*>& all, const int alpha, const int i) const
{
    double scatter = 0.;
    const double lam = lambda_[i];

    // scattering with other real pairs in this object
    for (int k=0; k< lambda_.size(); ++k) {
        // this excludes self-scattering of the real pair, absorbed into kinetic term
        if (k==i) continue;

        // note that if realpairs has a zero root, it's actually two roots very close to zero on either side.
        scatter += atan(lam-lambda_[k]) + atan(lam+lambda_[k]);
    }

    // scattering with others
    for (int beta=0; beta<all.size(); ++beta) {
        if (beta==alpha) continue;
        scatter += sum_re_atan(lam, all[beta]->hasOrigin(), all[beta]->getRealPairs(),
                           all[beta]->getImagPairs(), all[beta]->getComplexQuartets()).val();
    }

    // LHS == (N-1) atan (2x) == kinetic term + self-scattering
    return 0.5*tan( (0.5*PI*jx2_[i] + scatter) / (big_n_-1.)  ) ;
}



IterResult RealPairs::iterate(const vector<SymRoots*>& all, const int alpha,  double& convergence)
{
    vector<double> s (lambda_.size());
    //vector<double> lambda_next (lambda_.size());
    for (int i=0; i<lambda_.size(); ++i) {
        s[i] = step(all, alpha, i);
        z_next_[i] = secantStep(lambda_last_[i], lambda_[i], s_last_[i], s[i]);
        convergence += sq(z_next_[i] - lambda_[i]);
    }
    s_last_ = s;
    return IterResult::iter_ok;
}


void RealPairs::refresh()
{
    lambda_last_ = lambda_;
    lambda_ = z_next_;
}


IterResult RealPairs::initiate(const vector<SymRoots*>& all, const int alpha)
{
    for (int i=0; i<lambda_.size(); ++i) {
        if (jx2_[i] == 0)
            // initial value must not be zero
            // this will stay out of the way if the int*Piu or half-int*Pi for the others.
            lambda_[i] = 1e-2/big_n_;
        else
            lambda_[i] = tan(0.5*jx2_[i]*PI/big_n_);
    }
    lambda_last_ = lambda_;
    for (int i=0; i<lambda_.size(); ++i) {
        s_last_[i] = step(all, alpha, i);
    }
    lambda_ = s_last_;
    return IterResult::iter_ok;
}





/* NonsymRealRoots */




NonsymRealRoots::NonsymRealRoots(const int big_n, const vector<int>& jx2, const vector<double>& lambda)
    :   big_n_(big_n), jx2_(jx2), lambda_(lambda), lambda_last_(jx2.size()), s_last_(jx2.size()), z_next_(jx2.size())
{ }

NonsymRealRoots::NonsymRealRoots(const int big_n, const vector<int>& jx2)
      : big_n_(big_n), jx2_(jx2), lambda_(jx2.size()), lambda_last_(jx2.size()), s_last_(jx2.size()), z_next_(jx2.size())
{ }

int NonsymRealRoots::size() const
{   return lambda_.size();  }

vector<double> NonsymRealRoots::getRealRoots() const
{   return lambda_;  }

vector<double> NonsymRealRoots::getRapidities() const
{   return lambda_;  }

int NonsymRealRoots::stringLength() const
{   return 1; }


void NonsymRealRoots::add(const int jx2, const double lambda)
{
    jx2_.push_back(jx2);
    lambda_.push_back(lambda);
    lambda_last_.push_back(0.);
    s_last_.push_back(0.);
}


vector<CVal> NonsymRealRoots::getRoots() const
{
    // must copy to turn vector<double> into vector<complex<double>>
    vector<CVal> roots(lambda_.size());
    for (int i=0; i<roots.size(); ++i) {
        roots[i] = { lambda_[i], 0.,0,0. };
    }
    return roots;
}


double NonsymRealRoots::step(const vector<NonsymRoots*>& all, const int alpha, const int i) const
{
    double scatter = 0.;
    const double lam = lambda_[i];

    // scattering with other real pairs in this object
    for (int k=0; k< lambda_.size(); ++k) {
        // this excludes self-scattering of the real pair, absorbed into kinetic term
        if (k==i) continue;
        // note that if realpairs has a zero root, it's actually two roots very close to zero on either side.
        scatter += atan(lam-lambda_[k]);
    }

    // scattering with others
    for (int beta=0; beta<all.size(); ++beta) {
        if (beta==alpha) continue;
        const vector<double> other_x = all[beta]->getRealRoots();
        const vector<CVal> other_z = all[beta]->getComplexPairs();
        scatter += sum_re_atan(lam, other_x, other_z).val();
    }

    // LHS == (N-1) atan (2x) == kinetic term
    return 0.5*tan( (0.5*PI*jx2_[i] + scatter) / big_n_ ) ;
}



IterResult NonsymRealRoots::iterate(const vector<NonsymRoots*>& all, const int alpha, double& convergence)
{
    vector<double> s (lambda_.size());
    //vector<double> z_next (lambda_.size());
    for (int i=0; i<lambda_.size(); ++i) {
        s[i] = step(all, alpha, i);
        z_next_[i] = secantStep(lambda_last_[i], lambda_[i], s_last_[i], s[i]);
    }
    s_last_ = s;

    for (int i=0; i<lambda_.size(); ++i)
        convergence += sq(z_next_[i]-lambda_[i]);

    return IterResult::iter_ok;
}

void NonsymRealRoots::refresh()
{
    lambda_last_ = lambda_;
    lambda_ = z_next_;
}



IterResult NonsymRealRoots::initiate(const vector<NonsymRoots*>& all, const int alpha)
{

    for (int i=0; i<lambda_.size(); ++i) {
        // only initialise if this hasn't been done on construction
        if (lambda_[i]==0) {
            // nudge is essential if there's a zero q.n.
            if (jx2_[i]==0) lambda_[i] = + 1./big_n_;

            else lambda_[i] = 0.5*tan(0.5*jx2_[i]*PI/big_n_);
        }
    }
    lambda_last_ = lambda_;
    for (int i=0; i<lambda_.size(); ++i) {
        s_last_[i] = step(all, alpha, i);
    }

    lambda_ = s_last_;
    return IterResult::iter_ok;
}


