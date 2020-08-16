#include <iostream>
#include <vector>
using namespace std;

#include "scan.h"
#include "process.h"


int run(void)
{
	REAL delta, max_en;
	int number_sites, left_number_down, number_energy, function_index, base_element, number_holes;
	vector<int> base_vec (1, 0);
	cerr<<" function_index  delta  number_sites  left_number_down "<<endl;
	cin>> function_index >> delta >> number_sites >> left_number_down;
	cerr<< "input complete"<<endl;

	Chain* pdc = newChain(delta, number_sites);
	Base* pdb = newGroundBase(*pdc, left_number_down);
	State* p_dummy = newGroundState (*pdb);
	Quantity* p_quantity = newQuantity(function_index, *p_dummy);

	listContributionsPerMomentum (cout, *p_quantity, delta, number_sites, left_number_down);
}
