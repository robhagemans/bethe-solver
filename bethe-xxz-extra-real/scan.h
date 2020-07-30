#ifndef SCAN_H
#define SCAN_H

#include <vector>
#include <string>
// c-string library: for atoi(), atof() etc.
#include <string.h>
#include <sstream>
#include <fstream>

#include "generic.h"
#include "quantity.h"

using namespace std;

/** exceptions **/

extern const char* exc_InvalidIndex;


/** scan policies **/

// if set, solve inadmissible states too (requires symmetric solver)
#define ADMIT_INADMISSIBLE

// default threshold for string deviation
#define DEFAULT_DEVIATION_THRESHOLD 1e-20

/// TODO: do/do not reset-initial-values policy (global variable? create a class with static members?)



/** log files **/
// terminates log
#define CODA "------------------------------------------------------------------------------\n"

/** sentinel values **/
// don't set state id
#define NO_ID -1


/** add function class for scan function **/
// an add func takes form factor info and does something with it.

// empty add: do nothing
class AddFunc {
public:
	inline AddFunc(void) {};
	// do/record a calculation 
	virtual inline void operator() (Quantity& quantity) {};
	// the log
	virtual inline ostream& log(void) { return CNULL; };
	// the progress report
	virtual inline ostream& prog(void) { return CNULL; };
	// write to the log
	virtual inline void log (const string message) {};
	// write to the progress report
	virtual inline void prog (const string message) {};
	// log: declare state (don't finish line)
	virtual inline ostream& log (const State* const p_state) { return CNULL; };
	// log a message w/respect to a state
	virtual inline void log (const State* const p_state, const string message)	{};
	// log an exception w/respect to a state
 	virtual inline void log (const State* const p_state, const Exception exc) {};
};
	
	
// only log 
class OnlyLog: public AddFunc {	
protected:
	ostream& its_logstream;
	ostream& its_progstream;
public:	
	// don't log 
	inline OnlyLog (void): its_logstream(CNULL), its_progstream(CNULL) {};
	// write everything to one stream
	inline OnlyLog (ostream& logstream): its_logstream(logstream), its_progstream(logstream) {};
	// use a separate log and progress stream
	inline OnlyLog (ostream& logstream, ostream& progstream): its_logstream(logstream), its_progstream(progstream) {};
	// the log
	virtual inline ostream& log(void) { return its_logstream; };
	// the progress report
	virtual inline ostream& prog(void) { return its_progstream; };
	// write to the log
	virtual inline void log (const string message) { its_logstream << message <<endl; };
	// write to the progress report
	virtual inline void prog (const string message) { its_progstream << message <<endl; };
	// log: declare state (don't finish line)
 	virtual inline ostream& log (const State* const p_state) { return its_logstream<< name(p_state->p_base) <<"_state"<< p_state->id() <<": "; };
 	// log a message w/respect to a state
 	virtual inline void log (const State* const p_state, const string message)	{ log(p_state) << message <<endl; };
 	// log an exception w/respect to a state
 	virtual inline void log (const State* const p_state, const Exception exc) 	{ log(p_state) << exc <<endl; };
};


// gather the data in a matrix (binning)
// height: number of energy bins;   width: number of modes
class AddToMatrixFunc: public OnlyLog {
	Matrix<REAL>& table;
	REAL min_energy;
	REAL max_energy;
public:
	AddToMatrixFunc(
		Matrix<REAL>& table_form_factor, const REAL the_max_energy, const REAL the_min_energy=0.0, ostream& the_logstream=cerr)
		: OnlyLog(the_logstream), max_energy(the_max_energy), table(table_form_factor) {};
	virtual void operator() (Quantity& quantity);
};


// write to a form factor (S?ff) file
class AddToFile : public OnlyLog {
public:
	const static int output_precision;
protected:
	ostream& its_outfile; 
public:
	// don't log, just output
	AddToFile(void) : OnlyLog(), its_outfile(cout) { } ;
	// output all logging and progress report to single stream
	AddToFile(ostream& outfile, ostream& logstream) : OnlyLog(logstream), its_outfile(outfile) { its_outfile.precision(output_precision);  };
	// use a separate log and progress stream
	AddToFile(ostream& outfile, ostream& logstream, ostream& progress) : OnlyLog(logstream, progress), its_outfile(outfile) { its_outfile.precision(output_precision);  };
	
	virtual inline void operator() (Quantity& quantity)	{ 
		State* p_state = quantity.pRightState();
		its_outfile << name(p_state->p_base) <<SEP<< p_state->id() <<SEP
			<< quantity.mode() <<SEP<< quantity.momentum() <<SEP<< quantity.energy() <<SEP<< quantity.formFactor()  <<SEP
			<< p_state->convergence <<SEP<< p_state->stringDeviation() <<SEP<< p_state->devianceMagnitude() <<SEP<< p_state->iterations <<SEP<< p_state->newton_iterations  <<endl; 
	};
};


// read from a form factor (S?ff) file
inline istream& readFormFactor (
	istream& infile, 
	long long int& id, int& index_momentum, REAL& momentum, REAL& energy, REAL& form_factor, 
	REAL& convergence, int& iterations, int& newt_iter, REAL& deviation) 
	// read the base string into id (this will fail, but it's a dummy)
{	
	string dummy; REAL dummy_deviance_magnitude;
	return infile >> dummy >> id >>index_momentum >>momentum >>energy >>form_factor >>convergence >>deviation >> dummy_deviance_magnitude >>iterations >>newt_iter; 
}






/** scan result **/
struct ScanResult {
	int number_total;
	int number_calculated;
	int number_not_accepted;
	int number_forbidden;
	int number_equal;
	int number_deviated;
	int number_not_converged;
	int number_exceptions;
	REAL sum;
	ScanResult(void) : number_total(0), number_not_accepted(0),  number_forbidden(0), number_deviated (0), number_not_converged(0), number_equal(0),  number_calculated(0), number_exceptions(0), sum(0) {};
};



/** 'scan' one state **/
bool scanState(
	AddFunc& addfunc, 
	Quantity& quantity, 
	State* const p_state, 
	const REAL deviation_threshold,
	ScanResult& scan_result);



/** scan through a full base and calculate form factors **/
REAL scanBase (
	AddFunc& addFunc,
	Quantity& quantity,
	Base* p_base, 
	const REAL deviation_threshold = DEFAULT_DEVIATION_THRESHOLD,
	long long int id_start=NO_ID, long long int id_stop=NO_ID
);


#endif
