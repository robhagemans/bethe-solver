
/** class Engine **/

class Engine {
public:
	ScanIntervals* p_result;
	AddFunc& addfunc;
	Quantity& quantity;
	const Policy& policy;
	const AcceptFunc& accept;
	const REAL deviation_threshold;
public:
	Engine (AddFunc& the_addfunc, Quantity& the_quantity, const Policy& the_policy, const AcceptFunc& the_accept, const REAL the_deviation_threshold): addfunc(the_addfunc), quantity(the_quantity), policy(the_policy), accept(the_accept), deviation_threshold (the_deviation_threshold) 
	{	p_result = new ScanIntervals;	};
	~Engine (void) { delete p_result; };
	virtual void startScan (State* p_state, int start_hole, int stop_hole) =0;
	virtual ScanIntervals* collect (void) =0; 
	virtual void clear(void) =0;
};


class LocalEngine : public Engine {
public:
	Stopwatch& calculation_time;
public:
	LocalEngine (AddFunc& the_addfunc, Quantity& the_quantity, const Policy& the_policy, const AcceptFunc& the_accept, const REAL the_deviation_threshold, Stopwatch& the_calculation_time) : Engine (the_addfunc, the_quantity, the_policy, the_accept, the_deviation_threshold), calculation_time(the_calculation_time) {
		calculation_time.stop();
		calculation_time.reset();
	};
	inline void startScan (State* p_state, const int start_hole, const int stop_hole) { scanInnermostHole (p_state, start_hole, stop_hole); }
	inline ScanIntervals* collect (void) { return p_result; }; 
	inline void clear(void) { p_result->clear(); };
};


void scanRecursiveTEST (		
	Engine& engine,
	const ScanIntervals& all_scan_interval,	// all the ones done before
	const ScanIntervals& run_scan_interval,	
	ScanIntervals* p_last_scan_interval,	// expects the ones done the last run,
	REAL contribution_threshold,
	bool extend,
	Base* p_base);


	