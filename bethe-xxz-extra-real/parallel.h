#include "scan.h"
#include "pvm3.h"

/*
extern const int PVMEngine::msg_TaskExit = 100;
extern const int PVMEngine::msg_Result = 101;
extern const int PVMEngine::msg_Summary = 102;
*/
extern const char* exc_BadInput;
extern const char* exc_PVMSpawnFailed;
extern const char* exc_NoSummary;
extern const char* exc_Slave;
extern const char* exc_ChildNotInList;

class PVMEngine : public Engine {
public:
	vector<int> processes;
	ScanIntervals* p_result;
	string executable;
	static const int msg_TaskExit;
	static const int msg_Result;
	static const int msg_Summary;
	static const int msg_Exception;
public:
	PVMEngine (AddFunc& the_addfunc, const Quantity& the_quantity, const Policy& the_policy, const AcceptFunc& the_accept, const REAL the_deviation_threshold, string the_task) : Engine (the_addfunc, the_quantity, the_policy, the_accept, the_deviation_threshold), executable(the_task) { p_result = new ScanIntervals (&quantity.leftState()); };
	~PVMEngine (void) { pvm_exit(); };
	void startScan (State* p_state, int start_hole, int stop_hole);
	ScanIntervals* collect (void);
	void clear(void);
};



class PVMSendFunc : public AddFunc {
	int parent;
public:
	PVMSendFunc(void) : AddFunc() { if ((parent=pvm_parent()) == PvmNoParent) cerr<<"Warning: no pvm parent found"<<endl; }
	~PVMSendFunc(void) { pvm_exit(); };
	inline virtual void operator() (
		long long int id, int index_momentum, REAL momentum, REAL energy, REAL form_factor,
		int converged, int iterations, int newt_iter, REAL deviation);
/// these shouldn't be in this class...
	void send (Exception exc);
	void send (const Base& base, const Interval& summary);
};



/** test **/

bool scanTEST (
	AddFunc& addfunc,
	const Quantity& quantity,
	vector<Base*>& p_bases,
	const AcceptFunc& accept,
	const REAL deviation_threshold,
	const Policy policy,
	Stopwatch& calculation_time,
	const REAL contribution_threshold,
	const REAL threshold_factor);

