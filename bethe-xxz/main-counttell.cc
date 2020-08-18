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

	double delta = 1;
	int chain_length = 200;
	int number_down = 98;

	cout  << "delta  chain_length  number_down "<<endl;
	cin >> delta >> chain_length >> number_down;
	cout << "input complete"<<endl;

	Chain* p_chain = newChain (delta, chain_length, number_down);
	Base* p_ground_base = newGroundBase (*p_chain, number_down);

	vector<Base*> p_bases = allNewBases (p_chain, number_down, chain_length /* max string length*/, chain_length /*num particles*/, chain_length /*num spinons*/);

	long long int total_states = 0;
	for (int i=0;i<p_bases.size();++i) {
		long long int limid;
		cout<< name(*p_bases[i]) <<'\t'<< (limid = p_bases[i]->limId().back()) <<'\t'<< (total_states += limid) <<endl;
		cout<<p_bases[i]->limQuantumNumbers()<<endl;
		cout<<p_bases[i]->numberSlotsPerStringType()<<endl;
		cout<<p_bases[i]->limId()<<endl;
	}
	cout<<endl<< choose(chain_length, number_down) <<endl;
}


