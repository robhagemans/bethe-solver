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
#include "parallel.h"

int run(void)
{
	const char* here = "parallelslave";
	// read data from command line
	int start_hole, stop_hole;
	string chain_name, base_name, quantity_name, left_base_name;
	long long int id, left_id;
	string command_line = Command_Line.str();
	/// PVMSlaveReadArguments or some such?
	Command_Line >> start_hole >> stop_hole;
	///arguments<<*p_state<<SEP<<quantity;
	Command_Line >> chain_name >> base_name >> id >> quantity_name >> left_base_name >> left_id;
cerr<<id<<endl;
	PVMSendFunc pvm_send_func;
	try {
		// create everything
		Chain* p_chain = newChain(chain_name);
		Base* p_base = newBase (base_name, p_chain);
		Base* p_left_base = newBase (left_base_name, p_chain);
		if (!p_base || !p_left_base) throw Exception (here, exc_NullPointer, "base");
		State* p_state = newState (*p_base, id);
		State* p_left_state = newState (*p_left_base, left_id);
		if (!p_state || !p_left_state) throw Exception (here, exc_NullPointer, "state");
		p_left_state->solve();
		Quantity* p_quantity = newQuantity (quantity_name, *p_left_state);
		if (!p_quantity) throw Exception (here, exc_NullPointer, "quantity");
		Stopwatch calculation_time;
		/// policy, deviation_threshold and acceptfunc of calling function are currently ignored
		/// should really create PVMSlaveEngine or something, with sendfunc implied
		LocalEngine engine (pvm_send_func, *p_quantity, DEFAULT_POLICY, acceptAlways, DEFAULT_DEVIATION_THRESHOLD, calculation_time);
		// do the scan & send the results
		engine.startScan (p_state, start_hole, stop_hole);
		// send the interval record (there's one and only one)
		pvm_send_func.send (*p_base, engine.collect()->records[0].intervals[0]);
	}
	catch (Exception exc) {
		exc.remarks += "; parameters: " + command_line;
		pvm_send_func.send (exc);
	}

}


