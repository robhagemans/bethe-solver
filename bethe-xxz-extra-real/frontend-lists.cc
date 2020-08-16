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

int run(void)
{
	int number_sites;
	cerr<< "number sites"<<endl;
	cin>> number_sites;

	Chain* p_chain = new XXX_Chain (number_sites, 3);
	vector<int> base_vec (1, number_sites/4);
	Base* p_base1 = new XXX_Base (*p_chain, base_vec, 0, 0);
	Base* p_base2 = new XXX_Base (*p_chain, base_vec, 1, 0);
	Base* p_base3 = new XXX_Base (*p_chain, base_vec, 2, 0);
	--base_vec[0];
	Base* p_base4 = new XXX_Base (*p_chain, base_vec, 0, 1);
	Base* p_base5 = new XXX_Base (*p_chain, base_vec, 1, 1);
	Base* p_base6 = new XXX_Base (*p_chain, base_vec, 2, 1);
	cerr<< name(*p_base1) <<SEP<< p_base1->limId().back() <<endl;
	cerr<< name(*p_base2) <<SEP<< p_base2->limId().back() <<endl;
	cerr<< name(*p_base3) <<SEP<< p_base3->limId().back() <<endl;
	cerr<< name(*p_base4) <<SEP<< p_base4->limId().back() <<endl;
	cerr<< name(*p_base5) <<SEP<< p_base5->limId().back() <<endl;
	cerr<< name(*p_base6) <<SEP<< p_base6->limId().back() <<endl;
	listAccepted (0, p_chain, p_base1, number_sites/4,  0, 0, acceptHalfBrillouin);
	listAccepted (0, p_chain, p_base2, number_sites/4,  0, 0, acceptHalfBrillouin);
	listAccepted (0, p_chain, p_base3, number_sites/4,  0, 0, acceptHalfBrillouin);
	listAccepted (0, p_chain, p_base4, number_sites/4,  0, 0, acceptHalfBrillouin);
	listAccepted (0, p_chain, p_base5, number_sites/4,  0, 0, acceptHalfBrillouin);
	listAccepted (0, p_chain, p_base6, number_sites/4,  0, 0, acceptHalfBrillouin);
}
