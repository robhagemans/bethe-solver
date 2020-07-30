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
#include <csignal>

using namespace std;

#include "bethe.h"
#include "square.h"
#include "det.h"
#include "strip.h"
#include "exception.h"

#include "chain.h"
#include "base.h"
#include "state.h"


#include "recursive.h"



/** fatal signal: dump the intervals file and exit. **/

// globals for signal handler, is there a better way to do this?
Engine* SCANGINE = 0;
AddFunc* RECORDER = 0;
// avoid recursive call of signal handler 
volatile sig_atomic_t FATAL_ERROR_IN_PROGRESS = 0;
// see also glibc documentation, section 24.4.2
void fatalSignal(const int sig)
{
	// break any recursion
	if (FATAL_ERROR_IN_PROGRESS) raise (sig);
	FATAL_ERROR_IN_PROGRESS = 1;
	// dump intervals file.
	if (RECORDER) 
		RECORDER->prog()<< "caught signal "<< strsignal(sig) <<". writing intervals file."<<endl;
	if (SCANGINE) {
		SCANGINE->all_scan.merge(SCANGINE->last_scan);
		SCANGINE->all_scan.write(SCANGINE->quantity.fileName("int"));
	}
	if (RECORDER) {
		RECORDER->prog()<< "intervals file written."<<endl;
		RECORDER->prog()<< "calculation time "<<SCANGINE->calculation_time.humanReadable()<<endl;
		RECORDER->prog()<< "total time "<<SCANGINE->total_time.humanReadable()<<endl;
		RECORDER->prog()<< "good-bye."<<endl;
	}
	// reraise the signal with default handling, for correct exit
	signal (sig, SIG_DFL);
  	raise (sig);
}


int run(void)
{
	int quantity_type, chain_length, left_number_down, max_minutes;
	REAL delta, threshold;
	cout  << " quantity_type  delta  chain_length  left_number_down   threshold   max_minutes"<<endl;
	cin >> quantity_type >> delta >> chain_length >> left_number_down >> threshold >> max_minutes;
	cout << "input complete"<<endl;
	
	// create and solve the ground state
	State* p_ground_state = newGroundState (newGroundBase (newChain (delta, chain_length), left_number_down));
	p_ground_state->solve();
	// create the calculation quantity
	Quantity* p_quantity = newQuantity(quantity_type, p_ground_state);
	// prepare a list of bases
	vector<BaseData> bases = allNewBases (p_quantity, 3 /* max string length*/, 4 /*num particles*/, 4 /*num spinons*/);
	
	// output to .result file, errors to .log file, progress report to .prog or /dev/stderr
	ofstream outfile (p_quantity->fileName("result"), ios_base::app);
	ofstream logfile (p_quantity->fileName("log"), ios_base::app);
	ofstream progfile (p_quantity->fileName("prog"), ios_base::app);
	AddToFile add_file (outfile, logfile, progfile);
	
	add_file.prog()<<"recursive scan of form factor "<<p_quantity->name()<<" Delta="<<delta<<" N="<<chain_length<<" M_mu="<<left_number_down<<endl;
	add_file.prog()<<"maximum total time set to "<<max_minutes<<" min."<<endl<<endl;
	add_file.prog()<<"please fasten your seatbelts."<<endl;
	add_file.prog()<<"thank you for flying abacus airlines."<<endl<<endl;
	
	// solution policy parameters
	Policy policy = DEFAULT_POLICY;
	policy.max_iterations = 800;
	policy.precision = 1e-24;
	
	// setup all auxiliary stuff
	Engine scanner (add_file, *p_quantity, policy, /* deviation_threshold = */ 1e-10);
	scanner.max_seconds = 60 * max_minutes;
	
	// set pointers for signal handler
	SCANGINE = &scanner;
	RECORDER = &add_file;
	// register the signal handler
	signal(SIGTERM, fatalSignal); 
	signal(SIGINT, fatalSignal); 
	signal(SIGQUIT, fatalSignal); 
	signal(SIGABRT, fatalSignal); 
	signal(SIGSEGV, fatalSignal);
	
	// do the scan.
	scanner.scan (bases, threshold, 0.1);
	
	// sign off.
	add_file.prog() << "calculation time "<<scanner.calculation_time.humanReadable()<<endl;
	add_file.prog() << "total time "<<scanner.total_time.humanReadable()<<endl;
	for (int i=0;i<bases.size();++i) add_file.prog()<< name(bases[i])<<SEP<<scanner.all_scan.summary(bases[i])<<endl;	
	RECORDER->prog()<< "good-bye."<<endl;
	
	progfile.close();
	logfile.close();
	outfile.close();
}


