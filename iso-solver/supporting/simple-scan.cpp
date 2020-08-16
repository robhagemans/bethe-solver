#include <vector>
#include <string>
#include <iostream>

#include "exception.h"
#include "strip.h"
#include "configuration.h"
#include "iso-solver.h"

using namespace std;

const string exc_TooManyStrings = "impossible configuration - too many strings";

/** calculate the boundary for allowed values of the half-integer quantum number I **/
// for each string type j. lim_quantum == ceil(2*I_infty)
// so that quantum_number < lim_quantum (strictly lower!)
int limQuantumNumber (const int chain_length, const Configuration& config, const int type) 
{
	int sum_theta = 0;
	for (int k=0; k < config.numberTypes(); ++k) 
		// this sum_theta is Takahashi's formula, which has the correct counting:
		sum_theta += config.numberStringsOfType(k) * (1 + 2*(min(config.stringLength(type),config.stringLength(k))-1) + (type!=k) );
	
	return abs(chain_length - sum_theta);
}


int nextNumber (const int number)
{
    if (number>0) return -number;
    else return -number+2;
}

bool firstInType (const int chain_length, const Configuration& config, const int j, Strip<int>& quantum_number, 
             int max_alpha=-1 )
{
    if (max_alpha<0) max_alpha=config.numberStringsOfType(j);
    int next_qn = 1 - (config.numberStringsOfType(j)%2);
    for (int alpha=0; alpha < max_alpha; ++alpha) {
        // this happens only if it's a bad configuration - too many strings for a given type
        if (abs(next_qn) >= limQuantumNumber(chain_length, config, j)) return false;
    
        quantum_number.set(j, alpha, next_qn);
        next_qn = nextNumber(next_qn);
    }
    return true;
}

Strip<int> firstState (const int chain_length, const Configuration& config)
{
 	Strip<int> quantum_number (config);
 	
 	for (int j=0; j < config.numberTypes(); ++j)    
        if (!firstInType(chain_length, config, j, quantum_number)) {
            throw Exception("firstState", exc_TooManyStrings);
        }
 
    return quantum_number;
}

bool nextState (const int chain_length, const Configuration& config, Strip<int>& quantum_number)
{
 	for (int j=0; j < config.numberTypes(); ++j) {
 	    int lim_j = limQuantumNumber(chain_length, config, j);
        
        // set ground for types less than j
        for (int k=0; k<j; ++k) {
            firstInType(chain_length, config, k, quantum_number);
        }
               
        for (int alpha=0; alpha < config.numberStringsOfType(j); ++alpha) {
            int next = nextNumber(quantum_number(j, alpha));
            if (abs(next) >= lim_j)
	            continue;
	        else if (alpha<config.numberStringsOfType(j)-1 && next == quantum_number(j, alpha+1))
	            continue;
	        else {
	            quantum_number(j, alpha) = next;
	            // reset strings less than alpha
                firstInType(chain_length, config, j, quantum_number, alpha);
                
	            return true;
            } 
       
        }
    }
    return false;
}


bool goodConfiguration(const int chain_length, const Configuration& config) 
{
        
        // this code is a bit wasteful - must be a more efficient way...
        Strip<int> quantum_number (config);
        for (int j=0; j < config.numberTypes(); ++j)    
            if (!firstInType(chain_length, config, j, quantum_number)) return false;
        return true;    
}


vector<int> firstBase (const int number_down)
{
    return vector<int> (1, number_down);
}


bool nextBase (const int chain_length, vector<int>& base)
{
    vector<int> old_base = base;
    int number_down = 0;
    for (int j=0; j < base.size(); ++j) {
        number_down += base[j]*(j+1);
    }          
    for (int j=1; j < number_down; ++j) {
        if (j>=base.size()) base.resize(j+1);
        ++base[j];
        base[0] -= (j+1); // string length
                
        for (int k=1; k < j; ++k) { 
            base[0] += (k+1)*base[k];
            base[k] = 0;
        }
        if (base[0] >=0) { 
            if (goodConfiguration(chain_length, IsoConfiguration(base)))  return true; 
        }
        
        base = old_base;
    }      
    return false;
}




int main() 
{

    int chain_length = 48;
    int number_down = 4;
    cout<<"N M: ";
    cin >> chain_length >> number_down;
    cout<<"base: ";
    string base_str;
    vector<int> base;
    bool scan_all_bases = true;
    int base_size = 0;
    cin >> base_size;
    
    for (int i=0; i<base_size; ++i) {
        int next = 0;
        cin >> next;
        base.push_back(next);
        scan_all_bases = false;
    }
    
	// override standard policies and solve
	IsoBetheTakahashiSolver::precision = 1e-26*chain_length*chain_length;
	IsoBetheTakahashiSolver::max_iterations = 20000;
	
    cout.precision(8);
    cout<<showpos;
        
    try {
        int total_states = 0;
        int total_not_converged = 0;
        int base_nr = 0;
        if (scan_all_bases) base = firstBase(number_down);
        do {
            cout<<endl<<noshowpos<<"BASE#"<<base_nr<<": ("<<base<<")"<<endl;
            
            IsoConfiguration config (base);
            Strip<int> state = firstState(chain_length, config);    
            int num_states = 0;
            do {  
                cout<<noshowpos<<"#"<<base_nr<<"#"<<num_states<<showpos<<": 2I ["<<state<<"] ";
                ++num_states;
                
	            IsoSolver solver(chain_length, config, state);
	           	bool converged = solver.solve();
            	 
                
                cout<< " 2J ["<< solver.cleanQuantum2J() <<"] sums [";
                int sum = 0;
                for (int j=0; j<config.numberTypes(); ++j)
                for (int alpha=0; alpha<config.numberStringsOfType(j); ++alpha) {
                    sum += solver.sumOfBethe2xJ(j,alpha);
                    cout<<solver.sumOfBethe2xJ(j,alpha)<<";";
                }
                cout<<"] sum "<<sum; 
                
                cout<<" "<< solver.roots();	
                
            	if (!converged)  {
            	    cout<<" -- not converged";
            	    ++total_not_converged;
                }
                cout<<endl;
            } while (nextState(chain_length, config, state));  
            
            
            cout<< noshowpos<< num_states<< " states in base"<<endl;
            total_states += num_states;
            ++base_nr;
        } while (scan_all_bases && nextBase(chain_length, base));
 
        cout<<noshowpos;       
        cout << total_states<< " states scanned."<<endl;
        cout << total_not_converged<< " states failed."<<endl;
        
    }
    catch(Exception e) {
        cerr<<e<<endl;
    }  
    return 0;
}
