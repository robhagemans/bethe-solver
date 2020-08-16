#include <stdio.h>
#include <iostream>
#include <vector>
#include <complex>
#include <sstream>
#include <fstream>
#include <algorithm>
using namespace std;

//#include "scan.h"
#include "process.h"


struct EnergyFormFactor {
	REAL energy;
	REAL form_factor;
};
//
// determine if energy of a is less than that of b, for sorting
int AscendingEnergy (const EnergyFormFactor& a, const EnergyFormFactor& b)
{	return a.energy < b.energy; }



// generate structure factor file from form factor files that are found in the current directory
void makeSingleMomentumStructureFactorFileNoSmooth (const Quantity& quantity, REAL delta, int number_sites, int left_number_down, int given_index_momentum, REAL max_energy, int number_energy, REAL min_energy, REAL deviation_threshold)
{
	REAL sum = 0.0;
	REAL part_sum = 0.0;
	vector<REAL> structure_factor (number_energy);


	// construct the base of the file name
	stringstream description, infilename, outfilebase;
	description << nameFromValue(mon_Delta, delta) << nameFromValue(mon_Length, number_sites) << nameFromValue(mon_LeftDown, left_number_down, fieldWidth(number_sites));
	infilename << quantity.name() << "*" << description.str() << "*.result";
	outfilebase << quantity.sfName() << description.str() << "_k" << given_index_momentum<<"_nosmooth";
	// look for input files with these values
	vector<string> input_files = ls( infilename.str() );

	string logfile_name = uniqueName(outfilebase.str() + ".log");
	ofstream logfile (logfile_name.c_str());
	logfile <<"structure_factor " << quantity.sfName() <<endl;
	logfile <<"delta "<< delta <<endl;
	logfile <<"number_sites "<< number_sites <<endl;
	logfile <<"left_number_down "<< left_number_down <<endl;
	logfile <<"max_energy "<< max_energy <<endl;
	logfile <<"min_energy "<< min_energy <<endl;
	logfile <<"energy_resolution "<< number_energy <<endl;
	logfile <<"index_momentum "<<given_index_momentum<<endl;

	logfile.setf(ios_base::fixed);
	logfile.precision(FILE_PRECISION);
	cerr<<"creating vector"<<endl;
	vector<EnergyFormFactor> dataset;

	for (int i=0; i < input_files.size(); ++i) {

		ifstream infile(input_files[i].c_str());

		long long int id;
		int index_momentum, iter, newt;
		REAL energy, momentum, form_factor, convergence, deviation;
		while ( infile>>id>>index_momentum>>momentum>>energy>>form_factor>>convergence>>iter>>newt>>deviation) {
			if (index_momentum != given_index_momentum) continue;
			if (deviation > deviation_threshold) continue;

			EnergyFormFactor entry = { energy, form_factor };
			dataset.push_back(entry);
		}
		infile.close();
	}

	// sort the dataset in order of increasing energy
	sort (dataset.begin(), dataset.end(), AscendingEnergy);



	string outfile_name = uniqueName(outfilebase.str() + ".result");
	ofstream outfile (outfile_name.c_str());
	outfile.setf(ios_base::fixed);
	outfile.precision(FILE_PRECISION);

	// scale the data point by the distance to its neighbours
	REAL density = 1.0/ (dataset[1].energy - dataset[0].energy);
	outfile<< dataset[0].energy <<SEP<< 2.0*PI *dataset[0].form_factor * density <<SEP<<  2.0*PI*density/(1.0*number_sites) <<endl;
	if (dataset.size()) for (int i=1; i< dataset.size()-1; ++i) {
		REAL density = 2.0/ (dataset[i+1].energy - dataset[i-1].energy);
		outfile<< dataset[i].energy <<SEP<< 2.0*PI *dataset[i].form_factor * density <<SEP<< 2.0*PI*density/(1.0*number_sites) <<endl;
	}

	logfile << "total_sum "<<  sum/quantity.maxSum(number_sites, left_number_down) <<endl;
	logfile << "output_file " << outfile_name <<endl;
	logfile << CODA;
	outfile.close();
	logfile.close();
}


int run(void)
{
	REAL delta, max_en, min_en, peak_width;
	int number_sites, left_number_down, number_energy, function_index, mode;
	cout<<"quantity delta N M mode min_energy max_energy n_energy"<<endl;
	cin >> function_index >> delta >> number_sites >> left_number_down;
// 	cin >> mode;
	cin >> min_en >> max_en >> number_energy;
	cerr<< "input complete"<<endl;

	Chain* p_chain = newChain (delta, number_sites);
	Base* p_ground_base = newGroundBase (p_chain, left_number_down);
	State* dummy_state = newGroundState (p_ground_base);
	Quantity* quantity = newQuantity (function_index, dummy_state);
cerr<<"starting...";
	makeStructureFactorFile (*quantity, delta, number_sites, left_number_down, max_en, number_energy, min_en, 1e-20);
	return 0;
}
