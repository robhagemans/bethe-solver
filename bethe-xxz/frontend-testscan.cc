/*
 * properties of XXZ chains
 *
 * Rob Hagemans, instituut voor theoretische fysica, universiteit van Amsterdam
 *
 */

#include <vector>
#include <iostream>

using namespace std;

#include "exception.h"
#include "scan.h"

class AcceptBoundaryHolesFunc : public AcceptFunc {
public:
	int minNumberBoundaryHoles;
	int boundaryDistance;
	AcceptBoundaryHolesFunc (void);
	AcceptBoundaryHolesFunc (int minHoles, int distance);
	virtual bool operator() (State& state) const;
	virtual string name(void) const;
};

AcceptBoundaryHolesFunc::AcceptBoundaryHolesFunc(void) : minNumberBoundaryHoles(0), boundaryDistance(0)
{}

AcceptBoundaryHolesFunc::AcceptBoundaryHolesFunc(int minHoles, int distance) : minNumberBoundaryHoles(minHoles), boundaryDistance(distance)
{}

string AcceptBoundaryHolesFunc::name(void) const
{
	if (minNumberBoundaryHoles<=0) return "";
	stringstream sname;
	sname << "_boho-"<<minNumberBoundaryHoles<<"-"<<boundaryDistance;
	return sname.str();
}

bool AcceptBoundaryHolesFunc::operator() (State& state) const
{
	if (minNumberBoundaryHoles <= 0) return true;
	int number_boundary_particles =0;
	int number_boundary_holes =2*boundaryDistance; // one hole to the left, one to the right...
	/// assumes symmetric ground interval.
	int left_boundary_0 = state.base.numberStringsOfType(0);
	//int left_boundary_0 = state.base.numberDown();		// gives less speedup, higher contribution...

	for (int quantum=0; quantum < state.base.numberStringsOfType(0); quantum++) {

		int qn = state.quantum_number[0][quantum];
		// particles on the outer boundary
		if ((qn > (- left_boundary_0 - 2*boundaryDistance)) && (qn < -left_boundary_0)) ++number_boundary_particles; // outer left
		else if ((qn < ( left_boundary_0 + 2*boundaryDistance)) && (qn > left_boundary_0)) ++number_boundary_particles;	// outer right
		else if ((qn < (- left_boundary_0 + 2*boundaryDistance)) && (qn > -left_boundary_0)) --number_boundary_holes; // inner left
		else if ((qn > ( left_boundary_0 - 2*boundaryDistance)) && (qn < left_boundary_0)) --number_boundary_holes; // inner right
	}
	return (number_boundary_holes >= minNumberBoundaryHoles);
}



int run(void)
{
		int quantity_type;
		REAL delta;
		int chain_length;
		int number_types;
		int left_number_down;
		int number_particles, number_spinons, boundary;

		cout<< "  0 for Szff, nonzero for Smff >> ";
		cin>> quantity_type;

		cout<< "  anisotropy Delta >> ";
		cin>> delta;
		cout<< "  chain length N >> ";
		cin>> chain_length;

		cout<< "  number of down spins of left state M >> ";
		cin>> left_number_down;
		cout<< "  max number particles excluding spinons >> ";
		cin>> number_particles;
		cout<< "  max number spinons >> ";
		cin>> number_spinons;

//	cout<<"  boundary size >> ";
//	cin >> boundary;
boundary=0;

		Chain* p_chain = newChain (delta, chain_length, 4);
		Base* p_ground_base = newGroundBase (*p_chain, left_number_down);
		State* p_ground_state = newGroundState (*p_chain, *p_ground_base);
		Quantity* p_quantity = newQuantity(quantity_type, *p_ground_state);
		p_ground_state->solve();

		cout.precision(10);
		cout.setf(ios::fixed);
		REAL sum=0;

		vector< Base* > all_bases = allNewBases(p_chain, p_quantity->rightNumberDown(left_number_down), 2, number_particles, number_spinons+boundary); //2 is max string length
		for (int i=0; i< all_bases.size(); ++i)
			cout<< "s"<<all_bases[i]->numberSpinons() <<" h"<<all_bases[i]->numberHoles()
				<<" f"<<all_bases[i]->numberFreedoms()<<" "<< name(*(all_bases[i])) <<" id < "<<all_bases[i]->limId()<<endl;//.back()

		cout<<endl<<endl;
		int total_ids = 0;
		Stopwatch stopwatch;
		for (int i=0; i<all_bases.size(); ++i) {
			total_ids += all_bases[i]->limId().back()-1;
			//AcceptBoundaryHolesFunc acceptBoundaryHoles (all_bases[i]->numberSpinons() - number_spinons, int(ceil(0.25*(all_bases[i]->numberSpinons() - number_spinons))));
			REAL contrib = makeFormFactorFile (*p_quantity, p_chain, all_bases[i], left_number_down, DEFAULT_POLICY, 0,0, acceptAlways ); //acceptBoundaryHoles && acceptHalfBrillouin
			//REAL contrib = makeFormFactorFile (quantity, p_chain, all_bases[i], left_number_down);
			sum += contrib;
			double lap = stopwatch.lapUse();
			cout<< "cpl "<< contrib/lap <<" lap "<<lap<<" time "<<stopwatch.getUse()<<SEP<<SEP;
			cout<< "s"<<all_bases[i]->numberSpinons() <<" h"<<all_bases[i]->numberHoles() <<" f"<<all_bases[i]->numberFreedoms()<<" "<<name(*(all_bases[i])) ;
			cout<<" id < "<<all_bases[i]->limId().back()<<" total_ids "<< total_ids;
			cout<< SEP<<SEP<<SEP<< contrib <<SEP<< sum <<endl;
		}
		cout<<endl<<sum<<endl;

	return 0;
}
