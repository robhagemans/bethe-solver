/*
 * properties of XXZ chains
 *
 * Rob Hagemans, instituut voor theoretische fysica, universiteit van Amsterdam
 *
 */

#include <vector>
#include <string>
#include <iostream>
using namespace std;
#include "exception.h"
#include "scan.h"

int run(void)
{
	Policy policy = {1000, 1, 1e-28, 1e-2, 1e-2, 10.0, 1, 20}; // fewer newtons for large N
	vector<string> listfiles = ls("*.list");
	for (int i=0; i< listfiles.size(); ++i) {
		cerr<<"processing " <<listfiles[i]<<endl;
		makeFormFactorFileFromList (listfiles[i], policy);
	}
}


