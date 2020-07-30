#include <iostream>
#include <vector>
using namespace std;

#include "scan.h"
#include "process.h"


int run(void)
{
	REAL delta;
	int number_sites, left_number_down, function_index;
	cerr<<" function_index  delta  number_sites  left_number_down "<<endl;
	cin>> function_index >> delta >> number_sites >> left_number_down;
	cerr<< "input complete"<<endl;
	
	
	// TODO: have to separate Quantity from FormFactor or some such!!
	
	
	Chain* p_chain = newChain (delta, number_sites);
	Base* p_ground_base = newGroundBase (p_chain, left_number_down);		
	State* dummy_state = newGroundState (p_ground_base);
	Quantity* quantity = newQuantity (function_index, dummy_state);
	listContributions (cout, *quantity, delta, number_sites, left_number_down);
}
