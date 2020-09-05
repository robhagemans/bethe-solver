/* Rob Hagemans, instituut voor theoretische fysica, universiteit van Amsterdam */

#include <vector>
#include <iostream>
using namespace std;

#include "configuration.h"
#include "iso-solver.h"
//#include <initializer_list>

int main()
{
    IsoConfiguration config ({2,0,1});
	Strip<int> bt_numbers =  Strip<int> (config).set(0,0, 1).set(0,1, -5).set(2,0, 0);
	int number_sites = 10;


	IsoSolver solver(number_sites, config, bt_numbers);

	// override standard policies and solve
	IsoBetheTakahashiSolver::precision = 1e-26*number_sites*number_sites;
	IsoBetheTakahashiSolver::max_iterations = 20000;


	try {

    	if (!solver.solve())  cout<<"not converged"<<endl;

        cout.precision(15);

        cout<< solver.roots() <<endl;
        cout<< solver.cleanQuantum2J()<<endl;
        cout<<bt_numbers<<endl;

      	return 0;
    }
    catch (Exception e) {
        cout<<e<<endl;
        return 1;
    }
}
