#include "parallel.h"

/** PVMEngine **/


const int PVMEngine::msg_TaskExit = 100;
const int PVMEngine::msg_Result = 101;
const int PVMEngine::msg_Summary = 102;
const int PVMEngine::msg_Exception = 103;

const char* exc_BadInput = "malformed input";
const char* exc_PVMSpawnFailed = "failed to spawn PVM process";
const char* exc_Slave = "slave exception";
const char* exc_NoSummary = "slave did not send summary";
const char* exc_ChildNotInList = "killed child was not in process list (this shouldn't happen)";


void PVMEngine::startScan (
	State* p_state,
	int start_hole,
	int stop_hole)
{
	const char* here = "PVMEngine::startScan";
	if (!p_state) throw Exception (here, exc_NullPointer);
	int tid;
	stringstream arguments;
	arguments << start_hole <<SEP<< stop_hole <<SEP;
	///arguments<<*p_state<<SEP<<quantity;
	arguments << name(p_state->chain) +sep_Name <<SEP<< name(p_state->base)+sep_Name <<SEP<< p_state->id() <<SEP;
	arguments << quantity.name() <<SEP<< name (quantity.leftState().base)+sep_Name <<SEP<< quantity.leftState().id();
	// unload the stringstream into an array of c-strings
	const int child_argc =  8;
	char* child_argv[ child_argc+1 ];
	for (int i=0; i < child_argc; ++i) {
		string whats_read;
		if (!(arguments >> whats_read)) throw Exception (here, exc_BadInput);
		child_argv[i] = new char[whats_read.length()+1];
		strcpy (child_argv[i], whats_read.c_str());
	}
	child_argv[child_argc] = 0;
	// copy to a non-const c-string
	/// there must be a better way...
	char task[executable.length()];
	strcpy (task, executable.c_str());
	// spawn the task
	int numt = pvm_spawn(task, child_argv, PvmTaskDefault, (char*)0, 1, &tid);
	if (1 != numt) throw Exception(here, exc_PVMSpawnFailed);
	// ask to notify when exiting
	pvm_notify(PvmTaskExit, msg_TaskExit, 1, &tid);
	//pvm_notify(PvmTaskExit, msg_TaskExit, 1, &tid);
	processes.push_back(tid);
}

ScanIntervals* PVMEngine::collect(void)
{
	const char* here = "PVMEngine::collect";
	int  i=0;
	while (processes.size()) {
		// has our child died? then bury it
		// wait for the first exit message from any child
		int bufid = pvm_recv(-1, msg_TaskExit);
		if (bufid) {
			int tid = 0, bufid = 0;
			// get tid of dead child
			// (it is not the sender: a dying child sends no messages)
			int info = pvm_upkint(&tid, 1, 1);
			while (bufid = pvm_nrecv(tid, msg_Exception )) {
				int buffer_size=0, dummy;
				pvm_bufinfo(bufid, &buffer_size, (int*)0, (int*)0 );
				char buffer[buffer_size+1];
				// unpack the message to a c-string
				pvm_upkstr(buffer);
				// and log it
				///addfunc.log << buffer;
				throw Exception (here, exc_Slave, buffer);
			}
			// here we could use buf_info and initialise the buffer later...
			// receive all sent data from this process
			while (bufid = pvm_nrecv(tid, msg_Result )) {
				int buffer_size=0;
				pvm_bufinfo(bufid, &buffer_size, (int*)0, (int*)0 );
				char buffer[buffer_size+1];
				// unpack the message to a c-string
				pvm_upkstr(buffer);
				// convert this to something we can read
				stringstream converter;
				converter << buffer;
				// read the data
				///converter >> p_state;
				long long int id;
				int index_momentum, converged, iterations, newton_iterations;
				REAL momentum, energy, form_factor, deviation;
				readFormFactor (converter, id, index_momentum, momentum, energy, form_factor, converged, iterations, newton_iterations, deviation);
				// process the data
				///addfunc(*p_state);
				addfunc(id, index_momentum, momentum, energy, form_factor, converged, iterations, newton_iterations, deviation);
			}

			// receive the summary
			if (bufid = pvm_nrecv(tid, msg_Summary)) {
				int buffer_size=0;
				pvm_bufinfo(bufid, &buffer_size, (int*)0, (int*)0 );
				char buffer[buffer_size+1];
				// get the data
				pvm_upkstr(buffer);
				stringstream converter;
				converter << buffer;
				// read the data
				Interval done_interval;
				string base_name, chain_name;

				converter >> chain_name >> base_name >> done_interval;
				// insert the summary into our intervals list
				Chain* p_chain = newChain (chain_name);
				Base* p_base = newBase (base_name, p_chain);
				p_result->insert (p_base, done_interval);
				delete p_base;
				///TODO: copy chain as well in base copy thingy
				/// this still segfaults...
				///delete p_chain;
			}
			else throw Exception (here, exc_NoSummary);

			vector<int>::iterator index = find(processes.begin(), processes.end(), tid);
			if (index == processes.end()) throw Exception (here, exc_ChildNotInList); // this really shouldn't happen
			else processes.erase(index);
		}
		// try next child
	}
	return p_result;
}


void PVMEngine::clear(void)
{
	// kill all children, if any
	for (int i=0; i< processes.size(); ++i) pvm_kill(processes[i]);
	// clear process list
	processes.clear();
	// delete results
	p_result->clear();
}


/** PVMSendFunc **/


inline void PVMSendFunc::operator() (
	long long int id, int index_momentum, REAL momentum, REAL energy, REAL form_factor,
	int converged, int iterations, int newt_iter, REAL deviation)
{
	// format the output
	stringstream output;
	output.precision (FILE_PRECISION);
	output 	<< id <<SEP<< index_momentum<<SEP<< momentum <<SEP<< energy <<SEP<< form_factor <<SEP
			<< converged <<SEP<< iterations <<SEP<< newt_iter <<SEP<< deviation <<endl;

	if (parent == PvmNoParent) {
		// we're called locally. output to stdout.
		cout << output.str() <<endl;
	}
	else {
		// copy to avoid lack of const qualifier (wasn't there some cast for that?)
		char output_cstr [output.str().length()+1];
		strcpy (output_cstr, output.str().c_str());
		// clear send buffer
		pvm_initsend( PvmDataRaw );
		// pack the output
		pvm_pkstr(output_cstr);
		// and send to parent as result
		pvm_send(parent, PVMEngine::msg_Result);
	}
};


void PVMSendFunc::send (Exception exc)
{
	// send as one string
	string message = exc.location + ": " + exc.error + " " + exc.remarks;
	if (parent == PvmNoParent) {
		// we're called locally. output to stderr
		cerr << message <<endl;
	}
	else {
		char output_cstr [message.length()+1];
		strcpy (output_cstr, message.c_str());
		// clear send buffer
		pvm_initsend( PvmDataRaw );
		// pack the output
		pvm_pkstr(output_cstr);
		// and send to parent
		pvm_send(parent, PVMEngine::msg_Exception);
	}
};


void PVMSendFunc::send (const Base& base, const Interval& summary)
{
	// send as one string
	stringstream message;
	message << name(base.chain)+sep_Name <<SEP<< name(base)+sep_Name <<SEP << summary;
	if (parent == PvmNoParent) {
		// we're called locally. output to stderr
		cerr << message.str() <<endl;
	}
	else {
		char output_cstr [message.str().length()+1];
		strcpy (output_cstr, message.str().c_str());
		// clear send buffer
		pvm_initsend( PvmDataRaw );
		// pack the output
		pvm_pkstr(output_cstr);
		// and send to parent
		pvm_send(parent, PVMEngine::msg_Summary);
	}
};




/** should be moved back to scan.cc after PVMEngine is take out of this part **/

bool scanTEST (
	AddFunc& addfunc,
	const Quantity& quantity,
	vector<Base*>& p_bases,
	const AcceptFunc& accept,
	const REAL deviation_threshold,
	const Policy policy,
	Stopwatch& calculation_time,
	const REAL contribution_threshold,
	const REAL threshold_factor)
{
	const char* here = "scan";
	if (threshold_factor >= 1.0) throw Exception (here, exc_ThresholdOne);
	if (contribution_threshold >= 1.0) throw Exception (here, exc_ThresholdOne);

	// create empty interval object
	ScanIntervals all_scan_interval (&quantity.leftState());
	ScanIntervals run_scan_interval (&quantity.leftState());
	ScanIntervals* p_last_scan_interval = new ScanIntervals  (&quantity.leftState());
	ScanIntervals* p_scan_result = 0;

	//LocalEngine engine (addfunc, quantity, policy, accept, deviation_threshold, calculation_time);
	PVMEngine engine (addfunc, quantity, policy, accept, deviation_threshold, "parallelslave");

	addfunc.logstream <<"threshold 1.0"<<endl;

	for (int i=0; i< p_bases.size(); ++i) {
		addfunc.logstream <<"base "<<name(*p_bases[i])<<endl;
		scanRecursiveTEST(engine, all_scan_interval, run_scan_interval, p_last_scan_interval, 1.0, true, p_bases[i]);
	}
	all_scan_interval.merge(*engine.collect());
	engine.clear();
	addfunc.logstream <<"cumulative "<<all_scan_interval.numberCalculated()<<SEP<<all_scan_interval.contribution() <<endl;
	// steadily decrease threshold in steps of threshold_factor (<1.0)
	for (REAL threshold = threshold_factor; threshold >= contribution_threshold; threshold*=threshold_factor) {
		addfunc.logstream<< "threshold "<<threshold<<endl;
		sort (p_bases.begin(), p_bases.end(), BaseComparator(all_scan_interval));


		///
		for (int j=0; j<p_bases.size();++j) cerr<<name(*p_bases[j])<<SEP<<all_scan_interval.averageContribution(p_bases[j])<<endl;
		cerr<<endl;
		///

		// scan recursively over all bases given
		for (int i=0; i< p_bases.size(); ++i) {
			if (all_scan_interval.averageContribution(p_bases[i]) < threshold) continue;
			addfunc.logstream<< "base "<<name(*p_bases[i])<<endl;

			// extend scan
			scanRecursiveTEST (engine, all_scan_interval, run_scan_interval, p_last_scan_interval, threshold, true, p_bases[i]);
			p_scan_result = engine.collect();
			run_scan_interval.merge(*p_scan_result);

			do {
				// all_scan_interval should remain the same throughout one threshold-base-run
				// as it determines the way the tree is searched
				// last_scan_interval should always have only the results of the last sub-run
				// scan_interval should contain everything from this run.

				/// TODO: this solution is excessively ugly now that p_scan_result points to a member of engine
				delete p_last_scan_interval;
				p_last_scan_interval = p_scan_result;
				p_scan_result = new ScanIntervals  (&quantity.leftState());

				// sniff scan
				scanRecursiveTEST(engine, all_scan_interval, run_scan_interval, p_last_scan_interval, threshold, false, p_bases[i]);
				p_scan_result = engine.collect();
				run_scan_interval.merge(*p_scan_result);

			} while (p_scan_result->numberCalculated());

			all_scan_interval.merge(run_scan_interval);
			run_scan_interval.clear();
			engine.clear();
		}
		addfunc.logstream<< "cumulative "<<all_scan_interval.numberCalculated()<<SEP<<all_scan_interval.contribution() <<endl;
	}

	addfunc.logstream<< "threshold "<<contribution_threshold<<endl;
	for (int i=0; i< p_bases.size(); ++i) {
		addfunc.logstream<< "base "<<name(*p_bases[i])<<endl;

		scanRecursiveTEST(engine, all_scan_interval, run_scan_interval, p_last_scan_interval, contribution_threshold, true, p_bases[i]);
		p_scan_result = engine.collect();
		run_scan_interval.merge(*p_scan_result);

		do {
			// all_scan_interval should remain the same throughout one threshold-base-run
			// as it determines the way the tree is searched
			// last_scan_interval should always have only the results of the last sub-run
			// scan_interval should contain everything from this run.

			delete p_last_scan_interval;
			p_last_scan_interval = p_scan_result;
			ScanIntervals* p_scan_result = new ScanIntervals  (&quantity.leftState());

			scanRecursiveTEST(engine, all_scan_interval, run_scan_interval, p_last_scan_interval,  contribution_threshold, false, p_bases[i]);
			p_scan_result = engine.collect();
			run_scan_interval.merge(*p_scan_result);

		} while (p_scan_result->numberCalculated());

		all_scan_interval.merge(run_scan_interval);
		run_scan_interval.clear();
		engine.clear();
	}

//all_scan_interval.consolidate();
//cerr<<all_scan_interval<<endl;

	addfunc.logstream<< "total "<<all_scan_interval.numberCalculated()<<SEP<<all_scan_interval.contribution()<<endl;
}

