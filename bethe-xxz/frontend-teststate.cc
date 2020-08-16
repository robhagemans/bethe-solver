/*
 * properties of XXZ chains
 *
 * Rob Hagemans, instituut voor theoretische fysica, universiteit van Amsterdam
 *
 */

#include <math.h>
#include <complex>
#include <vector>
#include <iostream>

using namespace std;

#include "bethe.h"
#include "square.h"
#include "det.h"
#include "strip.h"
#include "exception.h"

#include "chain.h"
#include "base.h"
#include "state.h"


#include "scan.h"

int run(void)
{
		int quantity_type;
		cout<< "  0 for Szff, nonzero for Smff >> ";
		cin>> quantity_type;


		double zeta, delta;
		cout << "  anisotropy Delta >> ";
		cin >> delta;

		int chain_length;
		cout<< "  chain length N >> ";
		cin>>chain_length;

		Chain* p_chain = newChain (delta, chain_length);

		int state_id =0;
		int number_spinons=0;
		int number_types;
		int left_number_down;

		cout<<"  number of down spins of left state M >> ";
		cin >> left_number_down;

		cout<<"  number of non-real types in base >> ";
		cin>>number_types;
		vector<int>structure (p_chain->numberTypes());
		if (number_types+1> p_chain->numberTypes()) throw "Too many string types";

		int number_excited_rapidities=0;
		cout<<"  non-real part of base >> ";
		for (int i=1; i < number_types+1; ++i) {
			cin>>structure[i];
		}

		int number_infinite=0;
		cout<<"  number of infinite rapidities >> ";
		cin>>number_infinite;

		cout<<"  number of spinons  s >> ";
		cin >> number_spinons;


		cout<<endl<<endl;

		cout <<endl<<"(== for this chain ==)"<<endl;
		cout << "zeta " << p_chain->anisotropy <<endl;
		cout << "string lengths: ";
		for (int i=0; i < p_chain->numberTypes(); ++i) cout << p_chain->stringLength(i) << SEP;
		cout<<endl;

		Base* p_ground_base = newGroundBase (*p_chain, left_number_down);
		State* p_oregon = newGroundState (*p_ground_base);
		Quantity* p_quantity = newQuantity(quantity_type, *p_oregon);

		Base* p_base = newBase (*p_chain, p_quantity->rightNumberDown(left_number_down), structure, number_spinons, number_infinite);
		cout <<endl<<"(== for this base ==)"<<endl;

		string da_name = name(*p_chain, left_number_down) + name(*p_base);
		cout << "name: " << da_name <<endl;
		Chain* test_chain=0;
		Base* test_base = newBase(da_name, test_chain);
		cout << "check base name:" << name(*test_chain, left_number_down) + name(*test_base) <<endl;
		delete test_base;
		delete test_chain;

		cout << "lim quantum numbers: " << p_base-> limQuantumNumbers()<<endl;
		cout << "number of slots per string type: " << p_base->numberSlotsPerStringType () <<endl;
		cout << "number of slots per sector: " << p_base->numberSlotsPerSector () <<endl;
		cout << "number of down spins: "<< (p_base->numberDown())<<endl;
		cout << "number of spinons: "<< (p_base->numberSpinons())<<endl;
		cout << "number of holes: "<< (p_base->numberHoles())<<endl;
		cout << "minimum number holes: "<< p_base->minNumberHoles()<<endl;
		cout << "maximum number holes: "<< p_base->maxNumberHoles()<<endl<<endl;
		cout<<endl<<"  state number (0 <= id < "<<p_base->limId().back();
		cout<<")  id >> ";
		cin>> state_id;
		cout<<endl<<endl;

		cout<<"(== ground state ==)"<<endl;



		cout<< "quantum numbers " <<endl<< p_oregon->quantum_number << endl;
		cout <<"id " <<p_oregon->id()<<endl;
		if (!p_oregon->solve (DEFAULT_POLICY)) cout<<"NOT CONVERGED "<<endl;
		cout<< "rapidities "<<endl<<p_oregon->rapidity<<endl<<endl;
		cout<< "converged to "<< p_oregon->convergence << " after " << p_oregon->iterations <<" iterations and "<<p_oregon->newton_iterations <<" newton steps"<<endl;
		cout<< "momentum: " << p_oregon->momentum() <<endl;
		cout<< "energy: " << p_oregon->energy()<<endl;
		cout<< "norm: "<< p_oregon->norm() << endl;
		//cout<< "weight: "<< p_oregon->weight() <<endl;
		cout<< "deviation: " << p_oregon->stringDeviation() <<endl<<endl;
		cout<<endl;

		cout<<"(== state " << state_id << " ==)"<<endl;
		State* p_florida = newState (*p_base, state_id);

		cout << "quantum numbers " <<endl<< p_florida->quantum_number << endl;

		if (!p_florida->admissible()) cout<<"INADMISSIBLE STATE " <<endl;
		cout <<"recomputed id: " << p_florida->id()<<endl;

		//p_florida->solveBetheEquationInterpolateNewton (5000, iterations=0, 100, newtiter=0, 1e-28);
		Stopwatch stopwatch;
		if (!p_florida->solve (DEFAULT_POLICY)) cout<<"NOT CONVERGED "<<endl;
		cout << "converged to " << p_florida->convergence << " after " << p_florida->iterations <<" iterations and "<< p_florida->newton_iterations<<" newton steps."<<endl;
		cout << "iteration time "<<stopwatch.lapHumanReadable()<<endl;
		cout << "rapidities "<<endl<<p_florida->rapidity<<endl<<endl;
		cout << "momentum: " << p_florida->momentum() <<endl;
		cout << "energy: " << p_florida->energy()<<endl;
		cout << "norm: "<< p_florida->norm() << endl;

		if (p_quantity->name()==Transverse::ff_name)
				cout << "< GS | S- | "<<state_id<<" >^2 = " << p_florida->transverseFormFactor(*p_oregon) << endl;
		else 	cout << "< GS | Sz | "<<state_id<<" >^2 = " << p_florida->longitudinalFormFactor(*p_oregon) << endl;
		cout << "determinant time "<<stopwatch.lapHumanReadable()<<endl;


		cout<<endl;
		//cout << "weight: "<< p_florida->weight() <<endl;
		cout << "deviation: " << p_florida->stringDeviation() <<endl<<endl;
		cout << "total time "<<stopwatch.humanReadable()<<endl;
		cout<<endl;
}


