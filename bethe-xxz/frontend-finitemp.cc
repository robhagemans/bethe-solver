/*
 * Rob Hagemans, instituut voor theoretische fysica, universiteit van Amsterdam
 */

#include <vector>
#include <iostream>
using namespace std;

#include "exception.h"
#include "bethe.h"
#include "chain.h"
#include "base.h"
#include "state.h"
#include "scan.h"

#define MAX_STRING_LENGTH 2


int run(void)
{	
	
	int number_sites, number_energy;
	int function_index;
	int cutoff_states;
	REAL delta, max_energy, magnetic_field, beta;
	int from_number_holes, uptoinc_number_holes;
	
	cin >> function_index >> delta >> number_sites >> cutoff_states >> number_energy >> max_energy >> from_number_holes>> uptoinc_number_holes >> beta >> magnetic_field;
	cerr<<"input complete"<<endl;
	Chain* p_chain = newChain (delta, number_sites);	
	
	Matrix<REAL> szz (number_energy, p_chain->length());
	// we need all magnetisation subspaces for a finite-temp calculation
	for (int number_down_left=1; number_down_left <= number_sites/2; ++number_down_left) {
		
		/** find lowest states in this subspace **/
		vector<int> base_vec (1, number_down_left); // simplest base, only real roots.	
		/** watskeburt met de infinite_raps? **/
		
		vector< State* > lowest_states;
		// 1 or 2 holes must have lowest states.
		vector< Base* > all_bases = allNewBases (p_chain, number_down_left, MAX_STRING_LENGTH, 5, 0);
		
		for (int i=0; i<all_bases.size(); ++i) {
			for (long long int id = 0; id < all_bases[i]->limId().back(); ++id) {
				State* p_state;
				try {
					p_state = newState(*p_chain, *(all_bases[i]), id);
				} 
				catch (Exception exc) {
					if (exc.error == exc_Forbidden) continue; else throw;
				}
				try {
					if (!p_state->solve (DEFAULT_POLICY)) continue; // not converged? just try the next state.
				}
				catch (Exception exc) {
					cerr<<exc;
					continue;
				}
				REAL energy = p_state->energy();
				if (!lowest_states.size()) {
					lowest_states.push_back(p_state);
					continue;
				}
				if ((lowest_states.size() < cutoff_states ) || (energy < lowest_states.back()->energy())) {
					
					vector< State* >::iterator current_state = lowest_states.begin();
					
					while ((current_state != lowest_states.end()) && (energy >= (*current_state)->energy()) ) { 
						REAL current_energy = (*current_state)->energy();
						++current_state;
					}
					lowest_states.insert (current_state, p_state);	
					// shrink the vector to obey state number cutoff
					if (lowest_states.size() > cutoff_states) {

						delete lowest_states.back();
						lowest_states.pop_back();
					}
				}
				else {
					delete p_state; // weg errr--meeeee...
				}
			}
			// don't delete bases: they're being used for multiple states... just generate garbage.
		}
cerr<<"bases"<<endl;
for (int i=0; i<all_bases.size(); ++i) cerr<< *(all_bases[i]) <<endl;
cerr<<"lowest states: "<<endl;
for (int i=0; i< lowest_states.size(); ++i) 
	cerr<<i<<SEP<<(lowest_states[i])->base.numberHoles()
	<<" / "<<(lowest_states[i])->base<<SEP<<(lowest_states[i])->id()<<SEP<<(lowest_states[i])->indexMomentum()<<SEP<<(lowest_states[i])->energy()<<endl;

		/** scan through bases **/
		
		vector< Matrix<REAL>* > p_form_factor_matrices (cutoff_states);
		for (int i=0; i< lowest_states.size(); ++i) {
			p_form_factor_matrices[i] = new Matrix<REAL> (number_energy, p_chain->length());
			vector<Base*> all_right_bases = allNewBases (p_chain, rightNumberDown(function_index, number_down_left), MAX_STRING_LENGTH, uptoinc_number_holes, from_number_holes);
			for (int j=0; j < all_right_bases.size(); ++j) {
				addTableFormFactor (p_form_factor_matrices[i], function_index, p_chain, all_right_bases[j], number_down_left, max_energy, number_energy, lowest_states[i]);
				delete all_right_bases[j];
			}
		}
		
		
		for (int i=0; i < lowest_states.size(); ++i) {
			stringstream filename;
			filename << "test." << number_down_left << "." << i << ".result";
			ofstream outfile (filename.str().c_str());
			outfile.precision(10);
			outfile << (*p_form_factor_matrices[i]) << endl;
			outfile.close();				
			
			REAL left_energy = lowest_states[i]->energy() - number_down_left * magnetic_field; 
			(*p_form_factor_matrices[i]) *= 2.0 * PI * exp(- beta * left_energy);
			szz += (*p_form_factor_matrices[i])  ;
		}
		stringstream filename;
		ofstream outfile ("Szz-finitemp.result");
		outfile.precision(10);
		outfile<<szz<<endl;
		outfile.close();
		
	}
}


