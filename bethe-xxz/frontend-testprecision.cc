/*
 * properties of XXZ chains
 *
 * Rob Hagemans, instituut voor theoretische fysica, universiteit van Amsterdam
 *
 */

#include <math.h>
#include <complex>
#include <vector>
#include <iostream>

using namespace std;

#include "bethe.h"
#include "square.h"
#include "det.h"
#include "strip.h"
#include "exception.h"

#include "chain.h"
#include "base.h"
#include "state.h"


#include "scan.h"

int run(void)
{
		int quantity_type;
		cout<< "  0 for Szff, nonzero for Smff >> ";
		cin>> quantity_type;
		
		
		double zeta, delta;
		cout << "  anisotropy Delta >> ";
		cin >> delta;
		
		int chain_length;
		cout<< "  chain length N >> ";
		cin>>chain_length;

		Chain* p_chain = newChain (delta, chain_length);
		
		int state_id =0;
		int number_spinons=0;
		int number_types;
		int left_number_down;
		
		cout<<"  number of down spins of left state M >> ";
		cin >> left_number_down;
		
		cout<<"  number of non-real types in base >> ";
		cin>>number_types;
		vector<int>structure (p_chain->numberTypes());
		if (number_types+1> p_chain->numberTypes()) throw "Too many string types";
		
		int number_excited_rapidities=0;
		cout<<"  non-real part of base >> ";
		for (int i=1; i < number_types+1; ++i) cin>>structure[i];
		
		int number_infinite=0;
		cout<<"  number of infinite rapidities >> ";
		cin>>number_infinite;
		
		cout<<"  number of spinons  s >> ";
		cin >> number_spinons;
		
		Base* p_ground_base = newGroundBase (*p_chain, left_number_down);		
		State* p_oregon = newGroundState (*p_ground_base);
		Quantity* p_quantity = newQuantity(quantity_type, *p_oregon);
		Base* p_base = newBase (*p_chain, p_quantity->rightNumberDown(left_number_down), structure, number_spinons, number_infinite);		
		
		cout<<endl<<"  state number (0 <= id < "<<p_base->limId().back();
		cout<<")  id >> ";
		cin>> state_id;
		cout<<endl<<endl;
		cout.precision(20);
		
		State* p_florida = newState (*p_base, state_id);
		Policy policy = DEFAULT_POLICY;
		for (int i=0; i<50; ++i ) {
			policy.precision = pow (10, -1.0*i);
			p_oregon->setFreeRapidities();
			p_florida->setFreeRapidities();
			if (!p_oregon->solve (1000, policy.precision)) cout<<"GROUND STATE NOT CONVERGED "<<endl;
			if (!p_florida->admissible()) cout<<"INADMISSIBLE STATE " <<endl;
			if (!p_florida->solve (1000, policy.precision)) cout<<"NOT CONVERGED "<<endl;
			cout << policy.precision <<SEP<< p_oregon->convergence<<SEP<<p_florida->convergence <<SEP<< (*p_quantity) (*p_florida) <<endl;
		}
		
		

}


