#include <iostream>
using namespace std;

#include "exception.h"
#include "bethe.h"
#include "process.h"

int run(void)
{
	REAL delta, max_en;
	int number_sites, left_number_down, number_energy, function_index;
	cerr<<" function_index  delta  number_sites  left_number_down "<<endl;
	cin >> function_index >> delta >> number_sites >> left_number_down;
	cerr<< "input complete"<<endl;



	makeEqualTimeCorrelatorFile (function_index, delta, number_sites, left_number_down);

	return 0;
}
