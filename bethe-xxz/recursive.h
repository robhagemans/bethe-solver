#ifndef RECURSIVE_H
#define RECURSIVE_H

#include "generic.h"
#include "quantity.h"
#include "scan.h"


extern const char* exc_IndexTooHigh;
extern const char* exc_NegativeIndex;

extern const char* exc_NegativeInterval;
extern const char* exc_IntervalOverlap;

/** scan interval structure **/
class Interval {
public:
	long long int start;
	long long int stop;
	long long int number_calculated;
	long long int number_accepted;
	// interval's total, not average, contribution (i.e. sum divided by max_sum);
	REAL contribution;	
	Interval& Interval::operator+= (const Interval& plus);
	Interval& Interval::operator-= (const Interval& minus);
	Interval Interval::operator- (const Interval& minus) const;
};

extern const Interval NO_INTERVAL;

// outputting an Interval
ostream& operator<< (ostream& stream, const Interval& interval);

// inputting an Interval
istream& operator>> (istream& stream, Interval& interval) ;



class BaseRecord {
public:
	BaseData base_data;
	Interval summary;
	vector<Interval> intervals;
};

/** sort intervals according to begin value **/
inline bool intervalLessThan (const Interval& a, const Interval& b) 
{ 	return (a.start < b.start); }


/** class scanintervals **/

class ScanIntervals {
public:
	vector< BaseRecord > records;	
public:
	ScanIntervals() : records() {};
	~ScanIntervals(void) { clear(); };
	
	inline void clear(void) 
	{ 
// 		for (int i=0; i<records.size();++i)	delete records[i].p_base;
		records.clear(); 
	};
	
	void insert (const BaseData& base_data, Interval new_interval);
	void merge (const ScanIntervals& with_intervals);
	
	void consolidate (void);
	bool includes (const BaseData& base_data, const Interval& new_interval) const;
	bool includes (const BaseData& base_data, const long long int id) const;

	Interval summary(const BaseData& base_data, const long long int id_start=NO_ID, const long long int id_stop=NO_ID) const;
	
	
	REAL averageContribution(const BaseData& base_data, const long long int id_start=NO_ID, const long long int id_stop=NO_ID) const;
	inline long long int highestId (const BaseData& base_data, const long long int id_start=NO_ID, const long long int id_stop=NO_ID) const 
	{	long long int last = summary (base_data, id_start, id_stop).stop; if (NO_ID==last) return NO_ID; else return last-1; 	}; 
	
	inline REAL contribution(const BaseData& base_data, const long long int id_start=NO_ID, const long long int id_stop=NO_ID) const
	{	return summary (base_data, id_start, id_stop).contribution; };
	
	inline long long int numberCalculated (const BaseData& base_data, const long long int id_start=NO_ID, const long long int id_stop=NO_ID) const
	{	return summary (base_data, id_start, id_stop).number_calculated; };
	
	inline long long int numberAccepted (const BaseData& base_data, const long long int id_start=NO_ID, const long long int id_stop=NO_ID) const
	{	return summary (base_data, id_start, id_stop).number_accepted; };
	
	
	Interval summary(void) const;
	REAL contribution (void) const;
	long long int numberCalculated (void) const;
	long long int numberAccepted (void) const;
	inline REAL averageContribution(void) const
	{	return contribution()/(1.0*numberCalculated()); };
	
	// outputting a ScanIntervals
	ostream& write (ostream& stream) const; 
	void write (const string file_name) const;
	// inputting a ScanIntervals (merge with existing)
	istream& read (istream& stream);
	void read (const string file_name);

protected:	
	BaseRecord* ScanIntervals::recordForBase (const BaseData& base_data);
	const BaseRecord* ScanIntervals::recordForBase (const BaseData& base_data) const;
	BaseRecord* ScanIntervals::addBase (const BaseData& base_data);
};


/** sort bases **/
class BaseComparator {
public:
	ScanIntervals* intervals;
	BaseComparator (ScanIntervals* p_intervals) : intervals(p_intervals) {};
	inline int operator() (const BaseData& a, const BaseData&  b) const  {	
/*		
		if (intervals->numberCalculated(a)) {
			if (intervals->numberCalculated(b)) {
				// both have been checked; sort by highest contribution
				return intervals->averageContribution(a) > intervals->averageContribution(b); 
			}
			// only a has been checked; a is preferred	
			else return true;
		}
		// only b has been checked; b is preferred
		else if (intervals->numberCalculated(b)) return false;
		// none have been checked; prefer lowest number of particles
		else 
*/
		return a.numberFreedoms() < b.numberFreedoms();			
	};
};

void sortBases (vector<BaseData>& a, const int leftmost, const int rightmost, const BaseComparator& less_than);
	

/** class Engine **/

class Engine {
public:
	ScanIntervals last_scan;
	ScanIntervals all_scan;
	AddFunc& record;
	Quantity& quantity;
	const REAL deviation_threshold;
	Stopwatch calculation_time;
	Stopwatch total_time;
// running variables
//protected:
	int countdown;
	Base* p_base;
	State* p_state_start;
	// MAX_INT seconds is about 1600 years, so this should suffice
	int max_seconds; 	
public:
	// create engine
	Engine (AddFunc& the_record, Quantity& the_quantity, const REAL the_deviation_threshold);
	// destructor.
	~Engine(void) {};
	// scan recursively over a set of bases 
	bool scan (vector<BaseData>& p_bases, const REAL contribution_threshold, const REAL threshold_factor);
	
protected:
	// resursive part of scan algorithm 
	void scanRecursive  (const REAL contribution_threshold, int current_sector = NOT_SET, int current_particle = NOT_SET);
	// terminates recursion
	void scanInnermostHole (State* const p_state, const int start_hole, const int stop_hole);
	// scan a single state
	void scanState(State* const p_state, Interval& scan_result);
};




#endif
