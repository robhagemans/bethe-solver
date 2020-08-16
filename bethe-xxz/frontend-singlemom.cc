/*
 * properties of XXZ chains
 *
 * Rob Hagemans, instituut voor theoretische fysica, universiteit van Amsterdam
 *
 */

#include <vector>
#include <iostream>
using namespace std;
#include "exception.h"
#include "scan.h"

int run(void)
{
	int number_sites;
	int index_momentum;
	REAL momentum;
	int min_id, lim_id;
	cerr<< "number sites   momentum/PI  id>=  id<"<<endl;
	cin>> number_sites >> momentum >> min_id >> lim_id;
	index_momentum = int(floor((1.0*number_sites)*momentum/2.0));

	SingleMomentumFunc acceptMomentum (index_momentum);
	Chain* p_chain = new XXX_Chain (number_sites, 3);
	vector<int> base_vec (1, number_sites/4-1);
	Base* p_base1 = new XXX_Base (*p_chain, base_vec, 1, 0);
	Base* p_base2 = new XXX_Base (*p_chain, base_vec, 1, 1);
	Base* p_base3 = new XXX_Base (*p_chain, base_vec, 2, 1);
	makeFormFactorFile (0, p_chain, p_base1, number_sites/4, policy, min_id, lim_id, acceptMomentum);
	makeFormFactorFile (0, p_chain, p_base2, number_sites/4, policy, min_id, lim_id, acceptMomentum);
	makeFormFactorFile (0, p_chain, p_base3, number_sites/4, policy, min_id, lim_id, acceptMomentum);
	delete p_base1;
	delete p_base2;
	delete p_base3;
}
