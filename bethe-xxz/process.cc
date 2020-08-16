
#include "process.h"


/** process form factor files **/

// generate structure factor file from form factor files that are found in the current directory
// will work with both 'full' and 'hbz' files, as it simply ingnores the second half Brillouin zone and double counts the first.
void makeStructureFactorFile (const Quantity& quantity, REAL delta, int number_sites, int left_number_down, REAL max_energy, int number_energy, REAL deviation_threshold, REAL min_energy, const REAL peak_width)
{
	REAL sum = 0.0;
	REAL part_sum = 0.0;
	Matrix<REAL> structure_factor (number_energy, number_sites);

	// construct the base of the file name
	stringstream description;
	description << nameFromValue(mon_Delta, delta) << nameFromValue(mon_Length, number_sites) << nameFromValue(mon_LeftDown, left_number_down, fieldWidth(number_sites));

	// look for input files with these values
	vector<string> input_files = ls(quantity.name() + "*" + description.str() + "*.result");

	string logfile_name = uniqueName(quantity.sfName() + description.str() + ".log");
	ofstream logfile (logfile_name.c_str());
	logfile <<"structure_factor " << quantity.sfName() <<endl;
	logfile <<"delta "<< delta <<endl;
	logfile <<"number_sites "<< number_sites <<endl;
	logfile <<"left_number_down "<< left_number_down <<endl;
	logfile <<"max_energy "<< max_energy <<endl;
	logfile <<"min_energy "<< min_energy <<endl;
	logfile <<"energy_resolution "<< number_energy <<endl;
	logfile <<"peak_width_coefficient "<< peak_width <<endl;

	logfile.setf(ios_base::fixed);
	logfile.precision(FILE_PRECISION);

	// axis values
	ofstream xfile (uniqueName(quantity.sfName() + description.str() + ".xaxis").c_str());
	for (int i=0; i < number_sites; ++i) xfile << 2*PI*i/(1.0*number_sites) <<SEP;
	xfile <<endl;
	xfile.close();

	ofstream yfile (uniqueName(quantity.sfName() + description.str() + ".yaxis").c_str());
	for (int index_energy=0; index_energy < number_energy; ++index_energy) yfile << max_energy * index_energy / (1.0*number_energy) <<SEP;
	yfile<<endl;
	yfile.close();

	for (int i=0; i < input_files.size(); ++i) {
		REAL base_sum = 0.0;

		ifstream infile(input_files[i].c_str());
		long long int id;
		int index_momentum, iter, newt;
		REAL energy, momentum, form_factor, convergence, deviation;
		while ( readFormFactor(infile, id, index_momentum, momentum, energy, form_factor, convergence, iter, newt, deviation) ) {
			if (index_momentum > number_sites/2) continue; // if it isn't really the half brillouin zone, ignore other half.
			if (deviation > deviation_threshold) continue;

			// ignore everything outside the first half B.z.
			if ((index_momentum >=0) && (index_momentum <= number_sites/2)) {
				base_sum += form_factor;
				for (int index_energy=0; index_energy < number_energy; ++index_energy) {
					REAL probe_energy = min_energy + (max_energy - min_energy) * index_energy / (1.0*number_energy);
					structure_factor[index_energy][index_momentum] += 2.0*PI * form_factor * smoothDelta (probe_energy - energy, peak_width/(1.0*number_sites));
				}
			}
			// double count everything else except the B.z. boundary&center at k=0 and k=pi
			if ((index_momentum >0) && (index_momentum < number_sites/2)) {
				base_sum += form_factor;
				for (int index_energy=0; index_energy < number_energy; ++index_energy) {
					REAL probe_energy = min_energy + (max_energy - min_energy) * index_energy / (1.0*number_energy);
					structure_factor[index_energy][number_sites - index_momentum] += 2.0*PI * form_factor * smoothDelta (probe_energy - energy, peak_width/(1.0*number_sites));
				}
			}
		}
		infile.close();

		sum += base_sum;
		logfile << input_files[i] << " contribution "<< base_sum/quantity.maxSum(number_sites, left_number_down) << endl;
	}

	string outfile_name = uniqueName(quantity.sfName() + description.str() + ".result");
	ofstream outfile (outfile_name.c_str());

	outfile.setf(ios_base::fixed);
	outfile.precision(FILE_PRECISION);
	//outfile<<structure_factor<<endl;
	for (int i=0; i<number_energy; ++i){
		for (int j=0; j< number_sites; ++j) {
			outfile <<structure_factor[i][j]<<" ";
		}
		outfile<<endl;
	}

	logfile << "total_sum "<<  sum/quantity.maxSum(number_sites, left_number_down) <<endl;
	logfile << "output_file " << outfile_name <<endl;
	logfile << CODA;
	outfile.close();
	logfile.close();
}



// generate structure factor file from form factor files that are found in the current directory
void makeSingleMomentumStructureFactorFile (const Quantity& quantity, REAL delta, int number_sites, int left_number_down, int given_index_momentum, REAL max_energy, int number_energy, REAL min_energy, REAL deviation_threshold, const REAL peak_width)
{
	REAL sum = 0.0;
	REAL part_sum = 0.0;
	vector<REAL> structure_factor (number_energy);
	vector<REAL> density (number_energy);

	// construct the base of the file name
	stringstream description, infilename, outfilebase;
	description << nameFromValue(mon_Delta, delta) << nameFromValue(mon_Length, number_sites) << nameFromValue(mon_LeftDown, left_number_down, fieldWidth(number_sites));
	infilename << quantity.name() << "*" << description.str() << "*.result";
	outfilebase << quantity.sfName() << description.str() << "_k" << given_index_momentum;
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
	logfile <<"peak_width_coefficient "<< peak_width <<endl;
	logfile <<"index_momentum "<<given_index_momentum<<endl;

	logfile.setf(ios_base::fixed);
	logfile.precision(FILE_PRECISION);

	for (int i=0; i < input_files.size(); ++i) {
		REAL base_sum = 0.0;

		ifstream infile(input_files[i].c_str());
		long long int id;
		int index_momentum, iter, newt;
		REAL energy, momentum, form_factor, convergence, deviation;
		while ( readFormFactor(infile, id, index_momentum, momentum, energy, form_factor, convergence, iter, newt, deviation) ) {
			if (index_momentum != given_index_momentum) continue;
			base_sum += form_factor;
			for (int index_energy=0; index_energy < number_energy; ++index_energy) {
				REAL probe_energy = min_energy + (max_energy - min_energy) * index_energy / (1.0*number_energy);
				structure_factor[index_energy] += 2.0*PI * form_factor * smoothDelta (probe_energy - energy, peak_width/(1.0*number_sites));
				density[index_energy] += 2.0*PI * smoothDelta (probe_energy - energy, peak_width/(1.0*number_sites));
			}
		}
		infile.close();
		sum += base_sum;
		logfile << input_files[i] << " contribution "<< base_sum/quantity.maxSum(number_sites, left_number_down) << endl;
	}

	string outfile_name = uniqueName(outfilebase.str() + ".result");
	ofstream outfile (outfile_name.c_str());
	outfile.setf(ios_base::fixed);
	outfile.precision(FILE_PRECISION);
	for (int index_energy=0; index_energy < number_energy; ++index_energy)
		outfile<<min_energy + (max_energy - min_energy) * index_energy / (1.0*number_energy)<<SEP<<structure_factor[index_energy]<<SEP<<density[index_energy]/(1.0*number_sites)<<endl;

	logfile << "total_sum "<<  sum/quantity.maxSum(number_sites, left_number_down) <<endl;
	logfile << "output_file " << outfile_name <<endl;
	logfile << CODA;
	outfile.close();
	logfile.close();
}


// make a file of space-dependent equal-time correlators
// from form factor files found in the current directory
void makeEqualTimeCorrelatorFile (
	int function_index,
	REAL delta, int number_sites, int left_number_down,
	REAL deviation_threshold)
{
	REAL sum = 0.0;
	REAL part_sum = 0.0;
	vector<REAL> sum_form_factor (number_sites,0);

	// construct the base of the file name
	stringstream description;
	description << "*" << nameFromValue(mon_Delta, delta) << nameFromValue(mon_Length, number_sites) << nameFromValue(mon_LeftDown, left_number_down, fieldWidth(number_sites));

	string quantity_name = (function_index==0)?"Szff":"Smff";
	string sf_name = (function_index==0)?"Szz":"Smp";

	// look for input files with these values
	vector<string> input_files = ls(quantity_name + description.str() + "*.result");

	string logfile_name = uniqueName( sf_name + ("FT_t0_jdep_" + description.str()) + ".log");
	ofstream logfile (logfile_name.c_str());
	logfile <<"real_space_correlator " << sf_name <<endl;
	logfile <<"delta "<< delta <<endl;
	logfile <<"number_sites "<< number_sites <<endl;
	logfile <<"left_number_down "<< left_number_down <<endl;

	logfile.setf(ios_base::fixed);
	logfile.precision(FILE_PRECISION);

	for (int i=0; i < input_files.size(); ++i) {
		REAL base_sum = 0.0;
		bool half_zone = (string::npos != input_files[i].find ("hbz",0));
		ifstream infile(input_files[i].c_str());
		long long int id;
		int index_momentum, iter, newt;
		REAL energy, momentum, form_factor, convergence, deviation;
		while ( readFormFactor(infile, id, index_momentum, momentum, energy, form_factor, convergence, iter, newt, deviation) ) {
			if (abs(deviation) > deviation_threshold) continue;
			base_sum += form_factor;
			sum_form_factor[index_momentum] += form_factor;
			if ( half_zone && (index_momentum >0) && (index_momentum < number_sites/2)) {
				base_sum += form_factor;
				sum_form_factor[index_momentum] += form_factor;
			}
		}
		infile.close();

		sum += base_sum;
		logfile << input_files[i] << endl;//" contribution "//<< base_sum/quantity.maxSum(number_sites, left_number_down) << endl;
	}


	/* Fourier transform the sum vector */
	vector<REAL> real_space_correlator (number_sites, 0);
	for (int j=0; j < number_sites; ++j)
		for (int q=0; q < number_sites; ++q)
			real_space_correlator[j] += real( sum_form_factor[q] * exp(-I * (2.0*PI*q*j)/(1.0*number_sites) ))/(1.0*number_sites);

	string outfile_name = uniqueName(sf_name + "FT_t0_jdep_" + description.str() + ".result");
	ofstream outfile (outfile_name.c_str());

	outfile.setf(ios_base::fixed);
	outfile.precision(FILE_PRECISION);

	for (int j=0; j < number_sites; ++j) outfile << real_space_correlator[j] <<endl;
	outfile<<endl;
	outfile.close();

//	logfile << "total_sum "<<  sum/quantity.maxSum(number_sites, left_number_down) <<endl;
	logfile << "output_file " << outfile_name <<endl;
	logfile << CODA;
	logfile.close();
}


// make a file of time-dependent nth-neighbour correlators
// from form factor files found in the current directory
void makeNeighbourCorrelatorFile (const Quantity& quantity, const REAL delta, int number_sites, int left_number_down, REAL max_time, int number_time, int offset, REAL deviation_threshold)
{
	REAL sum = 0.0;
	REAL part_sum = 0.0;
	complex<REAL>* correlator = new complex<REAL> [number_time];

	// construct the base of the file name
	stringstream description;
	description << "*"<<nameFromValue(mon_Delta, delta) << nameFromValue(mon_Length, number_sites) << nameFromValue(mon_LeftDown, left_number_down, fieldWidth(number_sites));

	// look for input files with these values
	vector<string> input_files = ls(quantity.name() + description.str() + "*.result");

	stringstream fourier_name;
	fourier_name << "FT_j" << offset <<"_tdep_";

	string outfile_name = uniqueName(quantity.sfName() + fourier_name.str() + description.str() + ".result");
	ofstream outfile (outfile_name.c_str());

	outfile.setf(ios_base::fixed);
	outfile.precision(10);

	string logfile_name = uniqueName(quantity.sfName() + fourier_name.str() + description.str() + ".log");
	ofstream logfile (logfile_name.c_str());

	logfile <<"neighbour_correlator " << quantity.sfName() <<endl;
	logfile <<"neighbour "<< offset<<endl;
	logfile <<"delta "<< delta <<endl;
	logfile <<"number_sites "<< number_sites <<endl;
	logfile <<"left_number_down "<< left_number_down <<endl;
	logfile <<"max_time "<< max_time <<endl;
	logfile <<"time_resolution "<< number_time <<endl;

	logfile.setf(ios_base::fixed);
	logfile.precision(10);

	for (int i=0; i < input_files.size(); ++i) {
		REAL base_sum = 0.0;

		ifstream infile(input_files[i].c_str());
		long long int id;
		int index_momentum, iter, newt;
		REAL energy, momentum, form_factor, convergence, deviation;
		while (  readFormFactor(infile, id, index_momentum, momentum, energy, form_factor, convergence, iter, newt, deviation) ) {
			if (abs(deviation) > deviation_threshold) continue;

			base_sum += form_factor;
			for (int index_time=0; index_time < number_time; ++index_time) {
				REAL probe_time= max_time*index_time / (1.0*number_time);

				correlator[index_time] += form_factor * exp(-I * (probe_time*energy)) * exp(-I*(offset*momentum));
			}
		}
		infile.close();

		sum += base_sum;
		logfile << input_files[i] << " contribution "<< base_sum/quantity.maxSum(number_sites, left_number_down) << endl;
	}

	REAL normalise = sqrt(norm(correlator[0]));
	for (int index_time=0; index_time < number_time; ++index_time) {
		REAL probe_time= max_time*index_time / (1.0*number_time);
		correlator[index_time] /= normalise;
		outfile<< index_time<<SEP<<probe_time<<SEP<<real(correlator[index_time])<<SEP<<imag(correlator[index_time])<<SEP<<sqrt(norm(correlator[index_time]))<<endl;
	}
	outfile<<endl;
	outfile.close();

	logfile << "fraction_sum "<<  sum/quantity.maxSum(number_sites, left_number_down) <<endl;
	logfile << "output_file " << outfile_name <<endl;
	logfile << CODA;
	logfile.close();

	delete[] correlator;
}



// tabulate the contributions
// from form factor files found in the current directory
void listContributions (ostream& outfile, const Quantity& quantity, const REAL delta, const int number_sites, const int left_number_down)
{
	const char* here = "makeContribFile";
	//const int max_base_size=CUTOFF_TYPES;

	// construct the base of the file name
	stringstream description;
	description << "*" << nameFromValue(mon_Delta, delta) << nameFromValue(mon_Length, number_sites) << nameFromValue(mon_LeftDown, left_number_down, fieldWidth(number_sites));
	// look for input files with these values
	vector<string> input_files = ls(quantity.name() + description.str() + "*.result");

	outfile.precision(FILE_PRECISION);
	outfile.setf(ios_base::fixed);
	REAL sum=0.0;
	for (int i=0; i< input_files.size(); ++i) {
		ifstream infile (input_files[i].c_str());
/*
		// does the name contain "hbz"? if so, it is a half Brillouin zone.
		bool half_zone = (string::npos != input_files[i].find ("hbz",0));
*/
		bool half_zone= true; /// always assume half Brillouin zone
		long long int id;
		int index_momentum, iter, newt;
		REAL energy, momentum, form_factor, convergence, deviation;
		REAL base_sum = 0.0;
		while ( readFormFactor(infile, id, index_momentum, momentum, energy, form_factor, convergence, iter, newt, deviation) ) {
			/// note that assuming hbz on a non-hbz file will screw up.
			/// but I can't select for one half, since currently I get stuff on either half
			base_sum += form_factor;
			// if we only have a half zone, everything not on the zone boundary has to be counted double.
			if ( half_zone && (index_momentum != 0) && (index_momentum != number_sites/2)) base_sum += form_factor;
		}
		sum += base_sum;
		outfile<< base_sum/quantity.maxSum(number_sites, left_number_down) <<SEP;
		outfile<< sum /quantity.maxSum(number_sites, left_number_down)  <<SEP;
		outfile<< quantity.maxSum(number_sites, left_number_down)  <<SEP;
		outfile<< base_sum/quantity.maxSumHighestWeight(number_sites, left_number_down) <<SEP;
		outfile<< sum /quantity.maxSumHighestWeight(number_sites, left_number_down)  <<SEP;
		outfile<< quantity.maxSumHighestWeight(number_sites, left_number_down)  <<SEP;

		outfile<< input_files[i]<<SEP;
		if (half_zone) outfile << "(half Brillouin zone)";	// debug message
		outfile<<endl;
		infile.close();
	}
}


// make a momentum-contribution graph
// use only momenta in the lower half Brillouin zone if you plan to use hbz files as well.
void listContributionsPerMomentum (ostream& outfile, Quantity& quantity, const REAL delta, const int number_sites, const int left_number_down)
{
	const char* here = "makeContribFile";
	const int max_base_size=CUTOFF_TYPES;
	stringstream description;

	description << "*" <<nameFromValue(mon_Delta, delta) << nameFromValue(mon_Length, number_sites) << nameFromValue(mon_LeftDown, left_number_down, fieldWidth(number_sites) );

	vector<string> input_files = ls(quantity.name() + description.str() + "*result");
	outfile.precision(FILE_PRECISION);
	outfile.setf(ios_base::fixed);
	vector<REAL> sum (number_sites, 0.0);
	vector<REAL> first_moment (number_sites, 0.0);
	for (int i=0; i< input_files.size(); ++i) {
		ifstream infile (input_files[i].c_str());
		long long int id;
		int index_momentum, iter, newt;
		REAL energy, momentum, form_factor, convergence, deviation;
		while ( readFormFactor(infile, id, index_momentum, momentum, energy, form_factor, convergence, iter, newt, deviation) ) {
			sum[index_momentum] += form_factor;
			first_moment[index_momentum] += energy * form_factor;
		}
		//outfile<< base_sum/quantity.maxSum(number_sites, left_number_down) <<SEP;
		infile.close();
	}
	for (int q=0; q<number_sites; ++q)
	outfile
		<< q <<SEP<< sum[q] /quantity.maxSum(number_sites, left_number_down) <<SEP
		<< first_moment[q]/quantity.maxSum(number_sites, left_number_down) <<endl;
}


