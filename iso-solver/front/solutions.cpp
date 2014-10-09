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
/*
	// in c++98, use:
	//IsoConfiguration config = IsoConfiguration().set(0,1).set(2,1);
	
	// N=12 [ +8 ;; +4 ;]  [+6 -11 +6 -9] - this has apriori colliding bt numbers
	IsoConfiguration config ({1,0,1});
	Strip<int> bt_numbers =  Strip<int> (config).set(0,0, 8).set(2,0, 4);
	int number_sites = 12;
*/
/*	
	IsoConfiguration config = vector<int> ({2});
	Strip<int> bt_numbers =  Strip<int> (config).set(0,0, 7).set(0,1,5);
	int number_sites = 12;
	cout<<bt_numbers<<endl;	
*/

//	(0.688190960235473,0) (0.36327126400274,0)

/*	
	IsoConfiguration config = vector<int> ({0,1});
	Strip<int> bt_numbers =  Strip<int> (config).set(1,0, 6);
	int number_sites = 12;
	cout<<bt_numbers<<endl;	
*/	
	//(1.01846035050635,0.480831430911387) (1.01846035050635,-0.480831430911387)
/*
	IsoConfiguration config = vector<int> ({2,1});
	int number_sites = 12;
	Strip<int> bt_numbers =  Strip<int> (config);
	
	vector< complex<double> > roots = {0.688190960235473, 0.36327126400274, 1.01846035050635+0.480831430911387*I , 1.01846035050635-0.480831430911387*I};
	
    cout<<	IsoBetheTakahashiSolver::dirtyQuantumJforRoots(number_sites, roots);
    cout<<endl;
    return 0;
*/
   
	
	/*
	// this has a collapsing 4-string. I think there's a 3-string+1 solution nearby, but is it double counted?
	IsoConfiguration config ({1,0,1});
	Strip<int> bt_numbers =  Strip<int> (config).set(0,0, 0).set(2,0, 38);
	int number_sites = 48; //48
	*/
	
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
	
