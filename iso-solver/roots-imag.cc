
#include "roots.h"

using namespace std;



/** ImagPairs class */


ImagPairs::ImagPairs(const int big_n, const int dummy, const vector<int>& jx2)
    :   big_n_(big_n),
        level_(jx2.size()), im_lambda_(jx2.size()), im_lambda_last_(jx2.size()), s_last_(jx2.size()), z_next_(jx2.size())
{
    // with the exception of singular pairs,
    // the quantum number for the positive innermost imaginary root is:
    // - if no origin root, ie even number of particles,  (N-3)/2 and the perfect-string location is at 1.5i
    // - if origin root, ie odd number of particles, (N-2)/2 and the perfect-string location is at i.
    for (int i=0; i< jx2.size(); ++i) {
        level_[i] = 0.5*(big_n_-jx2[i]);
    }
}


ImagPairs::ImagPairs(const int big_n, const vector<double>& level)
    :   big_n_(big_n), level_(level), im_lambda_(level.size()), im_lambda_last_(level.size()), s_last_(level.size()), z_next_(level.size())
{   }

ImagPairs* ImagPairs::clone() const
{   return new ImagPairs(*this); }

int ImagPairs::size() const
{   return im_lambda_.size()*2; }


vector<IVal> ImagPairs::getImagPairs() const
{
    vector<IVal> pairs (im_lambda_.size());
    for (int i=0; i< im_lambda_.size(); ++i) {
        pairs[i] = { 0, im_lambda_[i] };
    }
    return pairs;
}


vector<CVal> ImagPairs::getRoots() const
{

    vector<CVal> roots (im_lambda_.size()*2);
    for (int i = 0; i< im_lambda_.size(); ++i) {
        roots[i] = { 0.,0.,0, im_lambda_[i] };
        roots[2*im_lambda_.size()-i-1] = { 0.,0.,0, -im_lambda_[i] };
    }
    return roots;
}




double ImagPairs::step(const vector<SymRoots*>& all, const int alpha, const int i) const
{
    const double lam = im_lambda_[i];
    double scatter = 0.;
    int real_part = 0;

    // scattering with other imag pairs in this object
    for (int j=0; j< im_lambda_.size(); ++j) {
        // this excludes self-scattering of the imag pair, absorbed into kinetic term
        if (i==j) continue;
        IVal l_m_lj = { 0, lam - im_lambda_[j] };
        IVal l_p_lj = { 0, lam + im_lambda_[j] };

        // this class doesn't keep track of pos&delta anyway, so no use splitting here
        scatter += 0.25*(im_atan_x4(l_m_lj) + im_atan_x4(l_p_lj));
        //real_part += re_atan_x4_div_pi(l_m_lj) + re_atan_x4_div_pi(l_p_lj);
        real_part += 2*(re_atan_x2_div_pi(l_m_lj) + re_atan_x2_div_pi(l_p_lj));

    }
    for (int beta=0; beta<all.size(); ++beta) {
        if (beta==alpha) continue;
        scatter += sum_im_atan({0, lam}, all[beta]->hasOrigin(), all[beta]->getRealPairs(),
                           all[beta]->getImagPairs(),all[beta]->getComplexQuartets(),
                           real_part);
    }
// real part should be zero to have a correct solution!
// if not, we should look at level-delta rather than level+delta

    // LHS == (N-1) atih (2y) == kinetic term + self scattering term.
    // this assumes wide pairs
    return 0.5/tanh( scatter / (big_n_-1.) ) ;
}


IterResult ImagPairs::iterate(const vector<SymRoots*>& all, const int alpha, double& convergence)
{
    IterResult result = IterResult::iter_ok;
    vector<double> s (im_lambda_.size());
    //vector<double> z_next (im_lambda_.size());

    for (int i=0; i<im_lambda_.size(); ++i) {
        s[i] = step(all, alpha, i);
        z_next_[i] = secantStep(im_lambda_last_[i], im_lambda_[i], s_last_[i], s[i]);

        // enforce minimum level from quantum number,
        // approach slowly when secant method is too enthusiastic.
        //double min_level = 0.5*(N_-jx2_[i]);
        if (z_next_[i] <= level_[i]) {
            z_next_[i] = 0.5* (level_[i] + im_lambda_[i]);
            // no convergence
            convergence += 1.;
            // signal that the result was obtained by constraining, so not valid
            result = IterResult::iter_constrained;
        }
        else
            convergence += norm(z_next_[i]-im_lambda_[i]);
    }

    im_lambda_last_ = im_lambda_;
    im_lambda_ = z_next_;

    s_last_ = s;
    return result;
}


void ImagPairs::refresh()
{
    //im_lambda_last_ = im_lambda_;
    //im_lambda_ = z_next_;
}

IterResult ImagPairs::initiate(const vector<SymRoots*>& all, const int alpha)
{
    // we need a small deviation that doesn't get above 1 if multiplied by number of pairs
    const double deviation = 1./big_n_;
    for (int i=0; i<im_lambda_.size(); ++i) {
        // with the exception of singular pairs,
        // the quantum number for the positive innermost imaginary root is:
        // - if no origin root, ie even number of particles,  (N-3)/2 and the perfect-string location is at 1.5i
        // - if origin root, ie odd number of particles, (N-2)/2 and the perfect-string location is at i.
        //double min_level = 0.5*(N_-jx2_[i]);

        // put the root at the string location with a small positive deviation.
        // note the distance between roots must never be exactly 1, things will break.
        im_lambda_[i] = level_[i]*(1.+2.*deviation);
        s_last_[i] = step(all, alpha, i);
        im_lambda_last_[i] = im_lambda_[i];
        im_lambda_[i] = level_[i]*(1.+deviation);
    }

    return IterResult::iter_ok;
}











/** ImagPairs2 class */


ImagPairs2::ImagPairs2(const int big_n, const vector<int>& pos)
    :   big_n_(big_n), pos_(pos), delta_(pos.size()), delta_next_(pos.size()), delta_last_(pos.size()), s_last_(pos.size())
{   }

int ImagPairs2::size() const
{   return pos_.size()*2; }

vector<IVal> ImagPairs2::getImagPairs() const
{
    vector<IVal> y (pos_.size());
    for (int i = 0; i< pos_.size(); ++i) {
        y[i] = { pos_[i] , delta_[i] };
    }
    return y;
}

vector<CVal> ImagPairs2::getRoots() const
{
    vector<CVal> root (pos_.size()*2);
    for (int i = 0; i< pos_.size(); ++i) {
        root[i] = {0.,0., pos_[i], delta_[i]};
        root[2*pos_.size()-i-1] = {0.,0., -pos_[i], -delta_[i]};
    }
    return root;
}


double ImagPairs2::step(const vector<SymRoots*>& all, const int alpha, const int i) const
{
    double scatter = 0.;
    int real_part = 0;
    IVal y = { pos_[i], delta_[i] };

    // scattering with other imag pairs in this object
    for (int j=0; j< pos_.size(); ++j) {
        IVal y_j = { pos_[j], delta_[j] };
        // this excludes self-scattering of the imag pair, absorbed into kinetic term
        if (i==j) continue;
        scatter += 0.25*(im_atan_x4(dif(y, y_j)) + im_atan_x4(sum(y, y_j)));
//        real_part += re_atan_x4_div_pi(dif(y, y_j)) + re_atan_x4_div_pi(sum(y, y_j));
          real_part += 2*(re_atan_x2_div_pi(dif(y, y_j)) + re_atan_x2_div_pi(sum(y, y_j)));

    }
    for (int beta=0; beta<all.size(); ++beta) {
        if (beta==alpha) continue;
        scatter += sum_im_atan( y, all[beta]->hasOrigin(), all[beta]->getRealPairs(),
                                all[beta]->getImagPairs(), all[beta]->getComplexQuartets(),
                                real_part);
    }

// real part should be zero to have a correct solution!
// if not, we should look at level-delta rather than level+delta
cerr<<real_part<<endl;

    // LHS == (N-1) atih (2y) == kinetic term + self scattering term.
    // this assumes wide pairs, but narrow pairs don't exist on the imag axis
    // (only the singular pair, which is exactly \pm 0.5)
    return -0.5*pos_[i] + 0.5/tanh( scatter / (big_n_-1.) ) ;
}


IterResult ImagPairs2::iterate(const vector<SymRoots*>& all, const int alpha, double& convergence)
{
    vector<double> s (pos_.size());
    for (int i=0; i<pos_.size(); ++i) {
        s[i] = step(all, alpha, i);
    }

    for (int i=0; i<pos_.size(); ++i) {
        delta_next_[i] = secantStep(delta_last_[i], delta_[i], s_last_[i], s[i]);

        // enforce minimum level from quantum number,
        // approach slowly when secant method is too enthusiastic.
        //double min_level = 0.5*(N_-jx2_[i]);
        /*if (z_next[i] <= level_[i]) {
            z_next[i] = 0.5* (level_[i] + im_lambda_[i]);
            // no convergence
            convergence += 1.;
        }
        else
        */
        convergence += norm(delta_next_[i]-delta_[i]);
    }

    s_last_ = s;
    return IterResult::iter_ok;
}


void ImagPairs2::refresh()
{
    delta_last_ = delta_;
    delta_ = delta_next_;
}


IterResult ImagPairs2::initiate(const vector<SymRoots*>& all, const int alpha)
{
    // we need a small deviation that doesn't get above 1 if multiplied by number of pairs
    const double deviation = 1./big_n_;
    for (int i=0; i<pos_.size(); ++i) {
        // with the exception of singular pairs,
        // the quantum number for the positive innermost imaginary root is:
        // - if no origin root, ie even number of particles,  (N-3)/2 and the perfect-string location is at 1.5i
        // - if origin root, ie odd number of particles, (N-2)/2 and the perfect-string location is at i.
        //double min_level = 0.5*(N_-jx2_[i]);

        // put the root at the string location with a small positive deviation.
        // note the distance between roots must never be exactly 1, things will break.
        delta_[i] = deviation*pos_[i];
        s_last_[i] = step(all, alpha, i);
        delta_last_[i] = delta_[i];
        delta_[i] = 0.5*deviation*pos_[i];
    }

    return IterResult::iter_ok;
}







/** Kite class **/

Kite::Kite (const int big_n, const int pos) : big_n_(big_n), pos_(pos)  { }


vector<CVal> Kite::getRoots() const
{   return {  {0.,0.,pos_,del_}, {0.,x_,0,0.}, {0.,0.,-pos_,-del_}, {0.,-x_,0,0.}  };  }

int Kite::size() const
{   return 4;   }

vector<double> Kite::getRealPairs() const
{   return { x_ };  }

vector<IVal> Kite::getImagPairs() const
{   return { {pos_, del_ } };  }


bool Kite::step(const vector<SymRoots*>& all, const int alpha, double& new_x, double& new_del) const
{
    double theta = 0.;
    double log_r = 0.;
    int re_atan_x4_div_pi = 0.;

    const IVal y = {pos_, del_};

    // scattering with others
    for (int beta=0; beta<all.size(); ++beta) {
        if (beta==alpha) continue;

        const bool other_o = all[beta]->hasOrigin();
        const vector<double> other_x = all[beta]->getRealPairs();
        const vector<IVal> other_y = all[beta]->getImagPairs();
        const vector<CVal> other_z = all[beta]->getComplexQuartets();

        // top root eq
        log_r += sum_im_atan(y, other_o, other_x, other_y, other_z, re_atan_x4_div_pi);
        // central eq
        theta += sum_re_atan(x_, other_o, other_x, other_y, other_z).val();
    }

cerr<<re_atan_x4_div_pi<<endl;
    if (re_atan_x4_div_pi%4) return false;

    // use top root equation to get norm
    // self scattering for top root eq, excluding both central +x and -x dangerous terms
    //log_r += 0.5*log( sq(x_) + sq(y.imag())  ); // 0.5*log(norm(z_+I);
    //log_r += zeta(x_, add_n_i(y, 1).imag()); // 0.5*log(norm(z_+I);
    log_r += 0.5*zeta_plus_x2( sum(x_, y) ); // 0.5*log(norm(z_+I);

    // top root kinetic term
    //log_r -= (big_n_-1) * 0.25*im_atan_x4(mul(y, 2));
    log_r -= (big_n_-1) * 0.25*im_atan_x4(mul(y, 2));


    // use central equation to get argument
    // self scattering for central eq, excluding both top and bottom dangerous terms
    //theta += xi(x_, add_n_i(y, 1).imag());
    theta += xi_plus( sum(x_, y) ).val();

    // central kinetic term
    theta -= (big_n_-1) * atan(2.*x_) ;

    new_x = exp(log_r)*cos(theta);
    new_del = exp(log_r)*sin(-theta);
    return true;
}


IterResult Kite::iterate(const vector<SymRoots*>& all, const int alpha, double& convergence)
{
    IterResult result = IterResult::iter_ok;

    double new_x = 0.;
    double new_del = 0.;
    if (!step(all, alpha, new_x, new_del))
        result = IterResult::iter_err_sign;

    // note that we do a secant multiplication with complex numbers here
    // which mixes up the x and y parts and seems to work better.
    // note that the adding and subtracting of pos could lead to a loss of precision
    complex<double> z (x_, 0.5*pos_+del_);
    complex<double> s (new_x, 0.5*pos_+new_del);

    z_next_ = secantStep(z_last_, z, s_last_, s);
    convergence += norm(z_next_-z);

    s_last_ = s;
    return result;
}


void Kite::refresh()
{
    z_last_ = x_ + I*(0.5*pos_+del_);
    x_ = real(z_next_);
    del_ = imag(z_next_) - 0.5*pos_;
}

void Kite::setInitial(const double x, const double del, const double x2, const double del2)
{
    z_last_ = x + I*(0.5*pos_+del);

    x_ = x2;
    del_ = del2;
}


IterResult Kite::initiate(const vector<SymRoots*>& all, const int alpha)
{
    double x1 = x_;
    double del1 = del_;
    if (z_last_==0.) {
        // first initial value not set, set now
        x_ = 0.01;
        del_ = 0.001;
    }
    else {
        x_ = real(z_last_);
        del_ = imag(z_last_) - 0.5*pos_;
    }

    double new_x = 0.;
    double new_del = 0.;
    step(all, alpha, new_x, new_del);
    s_last_ = new_x + I*(0.5*pos_+new_del);
    z_last_ = x_ + I*(0.5*pos_+del_);

    if(x1==0. && del1 ==0.) {
        // second inital value not set, use step value
        x_ = new_x;
        del_ = new_del;
    }
    else {
        x_ = x1;
        del_ = del1;
    }

    return IterResult::iter_ok;
}



