#include <cstdio>

#include <cmath>
#include "roots.h"

using namespace std;

bool not_equal(const double a, const double b) 
{
    return std::abs(a-b)>1e-14;
}


int main()
{
    //srand(0);
    HPIVal th = { 0, 0. };
    double th2 = 0.; 
    
    CVal z = {1e+100, 0., 0, 0.};
//        cout<<z.val()<<" "<<re_atan_x2(z).val()<<endl;
    th.add( re_atan_x2(z) );
    th2 += atan(z.real());

    cout<<th.n<<" "<<th.eps<<endl;
    cout<<0.5*th.val()<<" "<<sin(div(th,2))<<endl;
    cout<<th2<<" "<<sin(th2)<<endl;
    
    //th.add( re_atan_x2(z) );
    //th2 += atan(z.real());

    //th = mul(th, 1e+6);
    //th2 *= 1e+6;

    cout<<0.5*th.val()<<" "<<tan(div(th,2))<<endl;
    cout<<th2<<" "<<tan(th2)<<endl;

    for(int x=-10; x<=10; ++x) {
    
        for(int y=-4; y<=4; ++y) {
            CVal z = { 0.25*x, 0., 2*(y/2), 0.5*(y%2) };
            printf("%3i %3i %5.1f%+5.1fi at: %5.2f (%3i %5.2f) \n", x, y, z.real(), z.imag(), re_atan(z).real(), re_atan(z).n, re_atan(z).eps);
            
            printf(" xi+: %5.3f (%3i %+5.3f) ", xi_plus(z).real(), xi_plus(z).n, xi_plus(z).eps);
            printf(" xi-: %5.3f (%3i %+5.3f) \n", xi_minus(z).real(), xi_minus(z).n, xi_minus(z).eps);
            printf(" xi+: %5.3f (%3i %+5.3f) ", xi_plus_NEW(z).real(), xi_plus_NEW(z).n, xi_plus_NEW(z).eps);
            printf(" xi-: %5.3f (%3i %+5.3f) \n", xi_minus_NEW(z).real(), xi_minus_NEW(z).n, xi_minus_NEW(z).eps);

        }
    }

}
