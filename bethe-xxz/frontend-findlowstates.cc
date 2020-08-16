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


void findLowest (int function_index, REAL delta, int number_sites, int left_number_down, REAL max_energy, int number_energy)
{

	// construct the base of the file name
	stringstream description;
	description << name(delta, number_sites, left_number_down);

	// look for input files with these values
	vector<string> input_files = ls(Quantity_Name[function_index] + description.str() + "_h?_base*result");

	int lowest_id, lowest_index_momentum, lowest_iter, lowest_newt, lowest_converged;
	REAL lowest_energy=100000.0, lowest_momentum, lowest_form_factor;

	int second_id, second_index_momentum, second_iter, second_newt, second_converged;
	REAL second_energy=10000.0, second_momentum, second_form_factor;

	for (int i=0; i < input_files.size(); ++i) {
		REAL base_sum = 0.0;

		ifstream infile(input_files[i].c_str());
		int id, index_momentum, iter, newt, converged;
		REAL energy, momentum, form_factor;
		while ( infile >>id >>index_momentum >>momentum >>energy >>form_factor >>converged >>iter >>newt ) {
			if (!converged) continue;
			if (energy==0.0) continue;
			if ((index_momentum ==number_sites/2) && (energy <= lowest_energy)) {
				second_id = lowest_id;
				second_index_momentum = lowest_index_momentum;
				second_momentum =lowest_momentum;
				second_energy =lowest_energy;
				second_form_factor =lowest_form_factor;
				second_converged = lowest_converged;
				second_iter = lowest_iter;
				second_newt = lowest_newt;

				lowest_id = id;
				lowest_index_momentum = index_momentum;
				lowest_momentum =momentum;
				lowest_energy =energy;
				lowest_form_factor =form_factor;
				lowest_converged = converged;
				lowest_iter = iter;
				lowest_newt = newt;
			} else if ((index_momentum == number_sites/2-1) &&(energy<= second_energy)) {
				second_id = id;
				second_index_momentum = index_momentum;
				second_momentum =momentum;
				second_energy =energy;
				second_form_factor =form_factor;
				second_converged = converged;
				second_iter = iter;
				second_newt = newt;
			}

		}
		infile.close();


	}


	string outfile_name = uniqueName(Quantity_Name[function_index] + description.str() + ".lowest");
	ofstream outfile (outfile_name.c_str());

	outfile.setf(ios_base::fixed);
	outfile.precision(10);
	outfile	<< lowest_id<<SEP<< lowest_index_momentum <<SEP
				<< lowest_momentum <<SEP
				<< lowest_energy <<SEP
				<< lowest_form_factor <<SEP
				<< lowest_converged<<SEP
				<< lowest_iter <<SEP
				<< lowest_newt <<endl;
	outfile	<< second_id<<SEP<< second_index_momentum <<SEP
				<< second_momentum <<SEP
				<< second_energy <<SEP
				<< second_form_factor <<SEP
				<< second_converged<<SEP
				<< second_iter <<SEP
				<< second_newt <<endl;

	outfile.close();
}

int main(void)
{
	REAL delta, max_en;
	int number_sites, left_number_down, number_energy, function_index;
	cin >> function_index >> delta >> number_sites >> left_number_down;
	cin >> max_en >> number_energy;
	cerr<< "input complete"<<endl;

	findLowest (function_index, delta, number_sites, left_number_down, max_en, number_energy);

	return 0;
}
