#include <cstdlib>
#include <cmath>

#include "roots.h"

using namespace std;

bool not_equal(const double a, const double b) 
{
    return std::abs(a-b)>1e-14;
}

int main()
{
    srand(0);
    for (int i=0; i<1e6; ++i) {
        double a = 4*(0.5-rand()/double(RAND_MAX));
        double b = 4*(0.5-rand()/double(RAND_MAX));
        int c = 2-rand()/(RAND_MAX/4);
        double d = 4*(0.5-rand()/double(RAND_MAX));
        
        //cout<<a<<" "<<b<<" "<<c<<" "<<d<<endl;
        CVal z = {a,b,c,d};
        //cout<<z.val()<<endl;
        
        if (not_equal(re_atan_x2_CAUTION(z), re_atan_x2(z).val())) {
            cout<<"!"<<z.val()<<" "<<re_atan_x2_CAUTION(z)/PI<<" "<<re_atan_x2(z).val()/PI<<endl;
        }
        if (not_equal(xi_plus_CAUTION(z), xi_plus(z).val())) {
            cout<<"+" << z.val()<<" "<<xi_plus_CAUTION(z)<<" "<<xi_plus(z).val()<<endl; 
        }
        if (not_equal(xi_minus_CAUTION(z), xi_minus(z).val())) {
        
            cout<<"-" << z.val()<<" "<<xi_minus_CAUTION(z)<<" "<<xi_minus(z).val()<<endl;
        }
        
        //cout<<z.val()<<" "<<1+re_atan_x2(z)-re_atan_x2_NEW(z).val()<<" "<<re_atan_x2_NEW(z).n*0.5<<" "<<0.5*re_atan_x2_NEW(z).eps<<endl;
        //cout<<z.val()<<" "<<1+xi_plus(z)-xi_plus_NEW(z).val()<<" "<<xi_plus_NEW(z).n<<" "<<xi_plus_NEW(z).eps<<endl;
        //cout<<z.val()<<" "<<1+xi_minus(z)-xi_minus_NEW(z).val()<<" "<<xi_minus_NEW(z).n<<" "<<xi_minus_NEW(z).eps<<endl;
        
        //cout<<endl;        
    
    }    

}
