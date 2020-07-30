#include <stdio.h>
#include <iostream>
#include <vector>
#include <complex>
#include <sstream>
#include <fstream>
using namespace std;

#include "process.h"


// tabulate the contributions
// from form factor files found in the current directory
void makeDeltaE2 (ostream& outfile, const REAL delta, const int number_sites, const int left_number_down)
{
	const char* here = "calculateDeltaE2";
	const int max_base_size=CUTOFF_TYPES;
	
	// construct the base of the file name
	stringstream description;
	description << "*" << nameFromValue(mon_Delta, delta) << nameFromValue(mon_Length, number_sites) << nameFromValue(mon_LeftDown, left_number_down, fieldWidth(number_sites));
	// look for input files with these values
	vector<string> input_file_long = ls(Longitudinal::ff_name + description.str() + "*.result");
	vector<string> input_file_trans = ls(Transverse::ff_name + description.str() + "*.result");
	
	
	REAL sum = 0.0;
	for (int i=0; i< input_file_long.size(); ++i) {
		ifstream infile (input_file_long[i].c_str());
		// does the name contain "hbz"? if so, it is a half Brillouin zone.
		bool half_zone = (string::npos != input_file_long[i].find ("hbz",0));
		long long int id;
		int index_momentum, iter, newt;
		REAL energy, momentum, form_factor, deviation, convergence;
		REAL base_sum = 0.0;
		while ( readFormFactor(infile, id, index_momentum, momentum, energy, form_factor, convergence, iter, newt, deviation) ) {
			base_sum += sq(form_factor) * sq(cos(0.5*momentum)) ;
			// if we only have a half zone, everything not on the zone boundary has to be counted double.
			if ( half_zone && (index_momentum >0) && (index_momentum < number_sites/2)) 
				base_sum += sq(form_factor) * sq(cos(0.5*momentum));
		}
		sum += base_sum;
		infile.close();
	}
	
	for (int i=0; i< input_file_trans.size(); ++i) {
		ifstream infile (input_file_trans[i].c_str());
		// does the name contain "hbz"? if so, it is a half Brillouin zone.
		bool half_zone = (string::npos != input_file_trans[i].find ("hbz",0));
		long long int id;
		int index_momentum, iter, newt;
		REAL energy, momentum, form_factor, deviation, convergence;
		REAL base_sum = 0.0;
		while ( readFormFactor(infile, id, index_momentum, momentum, energy, form_factor, convergence, iter, newt, deviation) ) {
			base_sum += sq(form_factor) * sq(cos(0.5*momentum)) ;
			// if we only have a half zone, everything not on the zone boundary has to be counted double.
			if ( half_zone && (index_momentum >0) && (index_momentum < number_sites/2)) 
				base_sum += sq(form_factor) * sq(cos(0.5*momentum));
		}
		sum += 2.0 * base_sum;
		infile.close();
	}
	
	Chain* p_chain = newChain (delta, number_sites);
	REAL field = magneticField(*p_chain, left_number_down);
	
	outfile.precision(20);
	outfile.setf(ios::fixed);
	outfile << left_number_down <<SEP<< field <<SEP<< sum <<endl;
}



int run(void)
{
	REAL delta;
	int number_sites, left_number_down;
	cerr<<"delta N M"<<endl;
	cin >> delta >> number_sites >> left_number_down;
	cerr<< "input complete"<<endl;
	
	makeDeltaE2 (cout, delta, number_sites, left_number_down);
	
	return 0;
}
