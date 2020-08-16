#ifndef PROCESS_H
#define PROCESS_H

#include <string>
#include "scan.h"


#define DEFAULT_PEAK_WIDTH 0.8


/** output files **/
// precision() of output streams
#define FILE_PRECISION 15


/** structure factor files **/
void makeStructureFactorFile (
	const Quantity& quantity,
	REAL delta, int number_sites, int left_number_down,
	REAL max_energy, int number_energy, REAL min_energy=0.0,
	REAL deviation_threshold = BOOMBANOOMBA, const REAL peak_width = DEFAULT_PEAK_WIDTH);

void makeSingleMomentumStructureFactorFile (
	const Quantity& quantity,
	REAL delta, int number_sites, int left_number_down,
	int given_index_momentum,
	REAL max_energy, int number_energy, REAL min_energy=0.0,
	REAL deviation_threshold = BOOMBANOOMBA, const REAL peak_width = DEFAULT_PEAK_WIDTH);

/** structure factor fourier transforms **/
/*void makeEqualTimeCorrelatorFile (
	const Quantity& quantity,
	REAL delta, int number_sites, int left_number_down,
	REAL deviation_threshold = BOOMBANOOMBA);
*/
void makeEqualTimeCorrelatorFile (
	int function_index,
	REAL delta, int number_sites, int left_number_down,
	REAL deviation_threshold = BOOMBANOOMBA);

void makeNeighbourCorrelatorFile (
	const Quantity& quantity,
	const REAL delta, int number_sites, int left_number_down,
	REAL max_time, 	int number_time,
	REAL deviation_threshold = BOOMBANOOMBA);

/** contributions **/

// list the contributions to the sum rule of all files found
void listContributions (
	ostream& outfile,
	const Quantity& quantity,
	const REAL delta, const int number_sites, const int left_number_down);

void listContributionsPerMomentum (
	ostream& outfile,
	Quantity& quantity,
	const REAL delta, const int number_sites, const int left_number_down);

#endif
