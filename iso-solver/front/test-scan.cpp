#include <vector>
#include <iostream>

#include "basics.h"
#include "roots.h"

#include "takahashi-string.h"
#include "symmetric-string.h"
#include "nonsymmetric-string.h"

#include "solver.h"
#include "simplescan.h"

using namespace std;



struct StateInfo {
    vector<int> base;
    int state_n;
    vector<vector<int> > ix2;
    int iter;
    double conv;
    double betherr;
    IterResult result;
};


class Report {
public:
    static const int iterations=1000;
    static const double convergence;
    
    static int state_n_;
    static int base_n_;
    static vector< StateInfo > no_conv;
    static bool show_iter_;
        
    static bool stateReport(const int big_n, const int big_m, const vector<int>& base, const vector<vector<int> > ix2) 
    {
        ++state_n_;
        cout<<" base="<<base<<" ";
        cout<<showpos<<" 2I="<<ix2<<" "<<noshowpos;
        cout<<endl;
        return true;
    }


    static void base(const int big_n, const int big_m, const vector<int>& base, const int total_scanned, const int not_converged) 
    {
        state_n_=0;
        ++base_n_;
        cout<<endl;
        cout<<"---- summary for base "<<base<<" ----"<<endl ;
        cout<<total_scanned<<" states scanned"<<endl; 
        cout<<not_converged<<" not converged "<<endl; 
        cout<<"-------------------------------"<<endl ;
        cout<<endl;
    }

    static void final(const int big_n, const int big_m, const vector<int>& base, const int total_scanned, const int not_converged) 
    {
        
        cout<<noshowpos<<endl;
        cout<<"---- totals for this scan ----"<<endl;
        cout<<total_scanned<<" states scanned"<<endl; 
        cout<<not_converged<<" not solved: "<<endl; 
        for(int i=0; i<no_conv.size(); ++i) {
            cout<<" "<<no_conv[i].base<<"#"<<no_conv[i].state_n<<": "<<showpos<<no_conv[i].ix2<<noshowpos;
            cout<<" iter "<<no_conv[i].iter<<" conv "<<no_conv[i].conv<<" err "<<no_conv[i].betherr;
            cout<<"  "<<message(no_conv[i].result)<<endl;
        }
        cout<<noshowpos<<endl;
    }

    static void iteration(StateSolver* solver)
    {
        cout<< "iter "<< solver->iterations() <<" ";
        cout<<" conv "<< solver->convergence() <<" ";
        cout<<" rts "<< getValues(solver->getRoots()) <<" ";
        cout<<endl;
    }

    static bool stateSolve(const int big_n, const int big_m, const vector<int>& base, const vector<vector<int> > ix2) 
    {
        ++state_n_;
        cout<<endl<<noshowpos;
        cout<<" state "<<base<<"#"<<state_n_<<endl;
        cout<<showpos<<" 2I "<<ix2<<" "<<noshowpos<<endl;
        if (isSymmetric(ix2)) cout<<" symmetric"<<endl;    
            
        
        StateSolver* solver = show_iter_?
              newSolver (big_n, base, ix2, Report::iteration)
            : newSolver (big_n, base, ix2);
        IterResult converged = solver->solve(iterations, convergence);
        vector<CVal> roots = solver->getRoots();
        int iter = solver->iterations();
        double conv = solver->convergence();
        delete solver;
   
        cout.precision(15);  
        cout<< " iter " << iter <<endl;
        cout<< " conv " << conv <<endl;
        cout<< " roots " << getValues(roots) <<endl;
        cout.precision(4);
        
        const double betherr_limit = 1e-5;
        double betherr = 0;
        bool parityerr = false;
        cout<< " jx2 " << getJx2(big_n, roots, &betherr, &parityerr) <<endl;
        cout<< " err " << betherr <<endl;
        
        if (converged != IterResult::iter_ok) 
            cout<<"\033[1;31m "<<message(converged)<<" \033[0m\n"<<endl; //bold red text
        if (parityerr) 
            cout<<"\033[1;31m parity error \033[0m"<<endl; 
        if (betherr>betherr_limit) 
            cout<<"\033[1;31m large quantum number error \033[0m"<<endl;
        if (converged != IterResult::iter_ok || parityerr || betherr>betherr_limit ) {
           no_conv.push_back( { base, state_n_, ix2, iter, conv, betherr, converged } );
        }
        
        return (converged==IterResult::iter_ok);
    }

};

const double Report::convergence = 1e-20;
int Report::base_n_ = 0;
int Report::state_n_ = 0;
bool Report::show_iter_ = false;
vector< StateInfo > Report::no_conv;
        


int main() 
{

    int chain_length = 48;
    int number_down = 4;
    
    cout<<"N M: ";
    cin >> chain_length >> number_down;
    
    cout<<"base: ";
    //string base_str;
    vector<int> base;
    int base_size = 0;
    cin >> base_size;
    
    bool scan_all_bases = true;
    for (int i=0; i<base_size; ++i) {
        int next = 0;
        cin >> next;
        base.push_back(next);
        scan_all_bases = false;
    }
    
    int state=0;
    if (!scan_all_bases) {
        cout<<"state: ";   
        cin >> state;
    }
    
    SimpleScanner scanner(chain_length,  number_down);
 
    int no_conv=0;
    if (scan_all_bases) scanner.scan(Report::stateSolve, Report::base, Report::final);
    else if (!state) scanner.scanBase(Report::stateSolve, Report::base, no_conv, &base);
    else {
        Report::state_n_ = state-1;
        Report::show_iter_ = true;
        scanner.scanState(Report::stateSolve, &base, state);
    }
           
    return 0;
}
