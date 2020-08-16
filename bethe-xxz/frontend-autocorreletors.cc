#include <stdio.h>
#include <iostream>
#include <vector>
#include <complex>
#include <sstream>
#include <fstream>
using namespace std;

#include "bethe.h"
#include "square.h"
#include "det.h"
#include "strip.h"
#include "chain.h"
#include "base.h"
#include "state.h"

using namespace XXZ;
#include "scan.h"


#include "bethe.cc"
#include "square.cc"
#include "det.cc"
#include "strip.cc"
#include "chain.cc"
#include "base.cc"
#include "state.cc"
#include "scan.cc"



int main(void)
{
	REAL delta, max_time;
	int number_sites, left_number_down, number_time, function_index, deviation_threshold_power, offset;
	cin >> deviation_threshold_power;
	cin >> function_index >> delta >> number_sites >> left_number_down;
	cin >> max_time >> number_time >> offset;
	cerr<< "input complete"<<endl;

	makeNeighbourCorrelatorFile (function_index, delta, number_sites, left_number_down, max_time, number_time, offset,
		pow(10.0, -deviation_threshold_power));
}
