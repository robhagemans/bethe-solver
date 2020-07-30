/* Rob Hagemans, instituut voor theoretische fysica, universiteit van Amsterdam */

#include <vector>
#include <iostream>
using namespace std;

#include "generic.h"
#include "quantity.h"
#include "scan.h"






int run(void)
{	
	int max_number_sites;
	int number_down_left;
	
	Policy policy = {2000, 1, 1e-28, 1e-2, 1e-2, 10.0, 1, 20}; // fewer newtons for large N; DEFAULT_POLICY is better for small systems
	//Policy policy = {2000, 100, 1e-30, 1e-2, 1e-2, 10.0, 7, 20};

	cerr<<"N"<<endl;
	cin >> max_number_sites;
	cout <<"$N$ &"<<SEP<<"$\\lambda\\sub{num}$ &"<<SEP<<"$\\lambda_{N\\gg1}$ &"<<SEP<<"$\\delta\\sub{num}$ &"<<SEP<<SEP<<"$\\delta_{N\\gg1}$ &"<<"error \\\\ "<<endl;
	//<<SEP<<"<GS M=3|S+-|odsym>"<<endl;
	cout.precision(10);
	for (int number_sites=8; number_sites <= max_number_sites; number_sites +=2) {
		cout<<"$"<<number_sites<<"$ &"<<SEP;
		Chain* p_chain = newChain (1.0, number_sites, /* cutoff_types = */ 4);	
		Base* p_ground_base = newGroundBase (*p_chain, /*  number_down_left = */ 3);
		State* p_ground_state = newGroundState (*p_ground_base);
		
		
		p_ground_state->solve();
	
		//_R4_s0_base-1-0-1_id0:
		Base* p_base = newBase("_R4_s0_base-1-0-1", p_chain);
		State* p_state = newState(*p_base, 0);
		solve(p_state, policy,  1e-20);
		
		XXXDeviatedState* p_dev_state = (XXXDeviatedState*) p_state; // pointer cast. ouch.
		
		cout<<"$"<<abs(p_dev_state->rapidity(0,0))<<"$ &"<<SEP;
		cout<<"$"<<sqrt(12.0)*pow(3.0, -number_sites/2)<<"$ &"<<SEP;
		cout<<"$"<<abs(p_dev_state->deviance(1,0,0))<<"$ &"<<SEP;
		cout<<"$"<<24.0*(number_sites-1.0)*pow(3.0, -number_sites)<<"$ &"<<SEP;
		
		// this should be a method of XXX_State
		vector<complex<REAL> > bethe_j = p_dev_state->calculateBetheI();
		REAL error = 0.0;
		for (int i=0; i<bethe_j.size();++i) { 
			REAL two_j = 0.5*round(2.0*real(bethe_j[i]));
			error += norm(bethe_j[i] -  two_j);
// 			cout<<two_j<<SEP;
		}
		cout<<"$"<<sqrt(error)<<"$ \\\\"<<endl;
		
		Quantity* p_quantity = newQuantity(-1, *p_ground_state);
		p_quantity->setRightState(*p_state);
		try{
 			cout<<p_quantity->formFactor()<<SEP;
 		}
 		catch (...)
 		{
 			cout<<0<<SEP;
 		}
 		
 		p_dev_state->rapidity(0,0) = sqrt(12.0)*pow(3.0, -number_sites/2);
 		p_dev_state->rapidity(2,0) = -p_dev_state->rapidity(0,0);
 		p_dev_state->aberration(2,0,0) = p_dev_state->rapidity(0,0);
 		p_dev_state->aberration(2,0,2) = p_dev_state->rapidity(0,0);
 		p_dev_state->deviance(2,0,0) = 24.0*(number_sites-1.0)*pow(3.0, -number_sites);
 		p_dev_state->deviance(2,0,2) = -p_dev_state->deviance(2,0,0);
 		
 		// reset quantity
 		delete p_quantity;
 		p_quantity = newQuantity(-1, *p_ground_state);
 		p_quantity->setRightState(*p_state);
 		cout<<p_quantity->formFactor()<<endl;
 		
		delete p_state;
		delete p_base;
		delete p_quantity;
		delete p_ground_base;
		delete p_chain;
	}
	//return 0;

	
}
