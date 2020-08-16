#include "recursive.h"

const Interval NO_INTERVAL = {NO_ID, NO_ID, 0, 0, NOT_A_NUMBER };
const Interval ZERO_INTERVAL = {NO_ID, NO_ID, 0, 0, 0.0};

const char* exc_IndexTooHigh = "row or tableau index too high";
const char* exc_NegativeIndex = "negative row or tableau index";

const char* exc_NegativeInterval = "negative interval";
const char* exc_IntervalOverlap = "interval overlap";

// (approximate) maximum number of form factors per base calculated in initial phase
// TODO: move this somewhere else
#define MAX_INITIAL 4

/** output files **/
// precision() of output streams
#define FILE_PRECISION 15


/** interval addition **/
Interval& Interval::operator+= (const Interval& plus)
{
	number_calculated += plus.number_calculated;
	number_accepted += plus.number_accepted;
	contribution += plus.contribution;
	// start and stop are outer limits; not the full interval of summary has been calculated !
	if (NO_ID==start || plus.start < start) start = plus.start;
	if (NO_ID==stop || plus.stop > stop) stop = plus.stop;
	return (*this);
}

/** interval subtraction **/
// note that a +=b; a-=b is not the identity in this case (id numbers are left untouched...) !
Interval& Interval::operator-= (const Interval& minus)
{
	number_calculated -= minus.number_calculated;
	number_accepted -= minus.number_accepted;
	contribution -= minus.contribution;
	return (*this);
}

Interval Interval::operator- (const Interval& minus) const
{
	Interval copy = *this;
	return copy -= minus;
}


/** interval i/o **/
ostream& operator<< (ostream& stream, const Interval& interval)
{ 	return stream 	<< interval.start <<SEP<< interval.stop <<SEP<< interval.number_calculated <<SEP<< interval.number_accepted <<SEP<< interval.contribution; };

istream& operator>> (istream& stream, Interval& interval)
{	return stream 	>> interval.start >> interval.stop >> interval.number_calculated >> interval.number_accepted >> interval.contribution; };



/** find a base **/
BaseRecord* ScanIntervals::recordForBase (const BaseData& base_data) {
	const char* here = "ScanIntervals::recordForBase";
// 	if (!p_base) throw Exception (here, exc_NullPointer);
	int index_base = 0;
	for (; index_base < records.size(); ++index_base)
		if (records[index_base].base_data == base_data) break;
	if (records.size()==index_base) return 0;
	else return &(records[index_base]);
}

const BaseRecord* ScanIntervals::recordForBase (const BaseData& base_data) const {
	const char* here = "ScanIntervals::recordForBase() const";
// 	if (!p_base) throw Exception (here, exc_NullPointer);
	int index_base = 0;
	for (; index_base < records.size(); ++index_base)
		if (records[index_base].base_data == base_data) break;
	if (records.size()==index_base) return 0;
	else return &(records[index_base]);
}


/** add a base, return pointer to corresponding interval **/
BaseRecord* ScanIntervals::addBase (const BaseData& base_data)  {
	const char* here = "ScanIntervals::addBase";
// 	if (!p_base) throw Exception (here, exc_NullPointer);
	// copy base and chain so that they may be safely deleted
// 	Base* p_base_copy = p_base->clone();
// 	p_base_copy->p_chain = p_base->p_chain->clone();
	BaseRecord new_record = { BaseData(base_data), Interval(), vector<Interval>() };
	records.push_back (new_record);
	return &(records.back());
}

/** insert an interval **/
void ScanIntervals::insert (const BaseData& base_data, Interval new_interval)
{
	const char* here = "ScanIntervals::insert";
	// negative intervals are illegal
	if (new_interval.start > new_interval.stop) throw Exception (here, exc_NegativeInterval);
	// check for null base
// 	if (!p_base) throw Exception (here, exc_NullPointer);
	// empty intervals: don't bother
	if (new_interval.start == new_interval.stop) return;

	// find base
	BaseRecord* p_base_rec = recordForBase(base_data);
	// not found, add a copy of base to table
	// it is copied so that it may be simply deleted on destruction
	if (!p_base_rec) p_base_rec = addBase(base_data);

	// insert the interval at the back
	// should we do this in a more sorted way? or is it more efficient to just quicksort now and then?
	// pro always keeping the sort order: quicker to find whether interval is included
	// contra: lose time searching when inserting.
	p_base_rec->intervals.push_back(new_interval);
	// add to base summary
	p_base_rec->summary += new_interval;
};


/** check if an interval is included **/
bool ScanIntervals::includes (const BaseData& base_data, const Interval& new_interval) const
{
	const char* here = "ScanIntervals::includes";
	// empty intervals are always included
	if (new_interval.start == new_interval.stop) return true;
	// negative intervals are illegal
	if (new_interval.start > new_interval.stop) throw Exception (here, exc_NegativeInterval);
	// check

	// find base
	const BaseRecord* p_base_rec = recordForBase(base_data);
	if (!p_base_rec) return false;

	for (int i=0; i < p_base_rec->intervals.size(); ++i) {
		// NOTE: the intervals we add should be atomic, i.e. we assume no overlaps!
		// if this is not the case, there will be double counting.
		if ( (new_interval.start >= p_base_rec->intervals[i].start ) && (new_interval.start < p_base_rec->intervals[i].stop) ) {
			if ( (new_interval.stop > p_base_rec->intervals[i].start ) && (new_interval.stop <= p_base_rec->intervals[i].stop) ) return true;
			else throw Exception (here, exc_IntervalOverlap);
		}
	}
	return false;
}

/** check if an id is included **/
bool ScanIntervals::includes (const BaseData& base_data, long long int id) const
{
	const char* here = "ScanIntervals::includes";
// 	if (!p_base) throw Exception (here, exc_NullPointer);

	// find base
	const BaseRecord* p_base_rec = recordForBase(base_data);
	if (!p_base_rec) return false;
	// check
	for (int i=0; i < p_base_rec->intervals.size(); ++i)
		if ( (id >= p_base_rec->intervals[i].start ) && (id < p_base_rec->intervals[i].stop) )
			return true;
	return false;
}



/** merge two interval vectors & sort **/
void ScanIntervals::merge (const ScanIntervals& with_intervals)
{
	const char* here = "ScanIntervals::merge";
	// loop over bases in with_interval
	for (int i=0; i < with_intervals.records.size(); ++i) {
		// base in new vector to be inserted
		const BaseRecord* p_new_rec = &with_intervals.records[i];

		// find base in current interval vector
		BaseRecord* p_base_rec = recordForBase(p_new_rec->base_data);
		// not found, add base to table
		if (!p_base_rec) p_base_rec = addBase(p_new_rec->base_data);
		// insert all new elements at the back of this list
		// first resize to be big enough for all
		int old_size = p_base_rec->intervals.size();
		p_base_rec->intervals.resize(old_size + p_new_rec->intervals.size());
		// then copy the elements
		for (int k=old_size; k < p_base_rec->intervals.size(); ++k)
			p_base_rec->intervals[k] = p_new_rec->intervals[k-old_size];
		// and sort the result
		/// assuming the input is ordered, maybe it is quicker to insert one by one while looping the vector
		/// NOTE: currently we don't assume a sorted list, hence it's no use sorting it (unless consolidating)
		//sort (p_base_rec->intervals.begin(), p_base_rec->intervals.end(), intervalLessThan);
		// add summaries
		p_base_rec->summary += p_new_rec->summary;
	}
}


/** consolidate adjacent intervals **/
void ScanIntervals::consolidate (void)
{
	// loop over all bases
	for (int index_base=0; index_base < records.size(); ++index_base) {
		// get interval vector
		BaseRecord* p_base_rec = &(records[index_base]);
		// first quicksort the whole vector
		sort (p_base_rec->intervals.begin(), p_base_rec->intervals.end(), intervalLessThan);
		// then run through the vector finding adjacent intervals
		for (int index_interval=1; index_interval < p_base_rec->intervals.size(); ++index_interval) {
			// assume no overlaps. adjacent intervals?
			if ( p_base_rec->intervals[index_interval].start == p_base_rec->intervals[index_interval-1].stop )  {
				// merge the higher interval into the lower
				p_base_rec->intervals[index_interval-1] += p_base_rec->intervals[index_interval];
				// delete the higher
				p_base_rec->intervals.erase ( p_base_rec->intervals.begin() + index_interval );
				// and check at the same point again (WACHT tot het rode licht gedoofd is. Er kan NOG een trein komen.)
				--index_interval;

				/// NOTE: this erasing all the time may be time-inefficient
				/// if so, make a new vector, then copy (to reduce number of shifts)
				/// if necessary, a vector of pointers to vectors would be even faster, eliminating need to copy.
			}
		}
	}
}


Interval ScanIntervals::summary (const BaseData& base_data, const long long int id_start, const long long int id_stop) const
{
	const char* here = "ScanIntervals::summary";
// 	if (!p_base) throw Exception (here, exc_NullPointer);

	const BaseRecord* p_base_rec = recordForBase(base_data);
	if (!p_base_rec) return NO_INTERVAL;
	if ((NO_ID==id_start) && (NO_ID==id_stop)) return p_base_rec->summary;
	Interval the_summary = ZERO_INTERVAL;
	for (int i=0; i < p_base_rec->intervals.size(); ++i) {
		// do we find an interval that partly encompasses ours?
		// we include anything that overlaps fully or partly with the given [start, stop) interval
		if (	(NO_ID==id_start || p_base_rec->intervals[i].stop > id_start)
			 && (NO_ID==id_stop || p_base_rec->intervals[i].start < id_stop) ) {
			 the_summary += p_base_rec->intervals[i];
		}
	}
	return the_summary;
}


// calculate the average contribution of a given interval, NaN if not scanned
REAL ScanIntervals::averageContribution(const BaseData& base_data, const long long int id_start, const long long int id_stop) const
{
	Interval the_summary = summary (base_data, id_start, id_stop);
	return the_summary.contribution/(1.0*the_summary.number_calculated);
}


/** for a full base: **/

REAL ScanIntervals::contribution (void) const{
	REAL contrib = 0.0;
	for (int i=0; i<records.size(); ++i)
		contrib += records[i].summary.contribution;
	return contrib;
}

long long int ScanIntervals::numberCalculated (void) const {
	long long int number_calculated = 0;
	for (int i=0; i<records.size(); ++i)
		number_calculated += records[i].summary.number_calculated;
	return number_calculated;
}

long long int ScanIntervals::numberAccepted (void) const {
	long long int number_accepted = 0;
	for (int i=0; i<records.size(); ++i)
		number_accepted += records[i].summary.number_accepted;
	return number_accepted;
}

Interval ScanIntervals::summary (void) const {
	Interval result = ZERO_INTERVAL;
	for (int i=0; i<records.size(); ++i)
		result += records[i].summary;
	return result;
}


// outputting a ScanIntervals
ostream& ScanIntervals::write (ostream& stream) const
{
	stream.precision(FILE_PRECISION);
	for (int r=0; r < records.size(); ++r)
		for (int i=0; i < records[r].intervals.size(); ++i)
			stream << name(records[r].base_data) <<SEP << records[r].intervals[i] <<endl;
	return stream;
}
void ScanIntervals::write (const string file_name) const
{
	ofstream file_out (file_name.c_str());
	write (file_out);
	file_out.close();
}

// inputting a ScanIntervals
istream& ScanIntervals::read (istream& stream)
{
	string base_name;
	Interval interval;
	while (stream >> base_name >> interval) {
		BaseData base_data = readBaseData (base_name);
		insert (base_data, interval);
	}
	return stream;
}
void ScanIntervals::read (const string file_name)
{
	ifstream file_in (file_name.c_str());
	read (file_in);
	file_in.close();
}




/** range-checked quicksort (JS's version+ range check) **/
// std::sort seems to access element -1 (and so does this routine if we don't range check)
void sortBases (vector<BaseData>& a, const int leftmost, const int rightmost, const BaseComparator& less_than)
{
	static BaseData pivot;
	static int j;
	int i;
	const int array_size = a.size();

	if (rightmost > leftmost) {
		pivot = a[rightmost]; i = leftmost-1; j = rightmost;
		for (;;) {
			// NOTE: range check!
			while (++i<array_size &&  less_than(a[i], pivot));
			// NOTE: avoid access of element -1 !!
			while (--j>=0 && less_than(pivot, a[j]));
			if (i >= j) break;
			std::swap(a[i], a[j]);
		}
		std::swap(a[i], a[rightmost]);
		sortBases(a, leftmost, i-1, less_than);
		sortBases(a, i+1, rightmost, less_than);
	}
}




/** class Engine **/

/** constructor **/
Engine::Engine (AddFunc& the_record, Quantity& the_quantity,  const REAL the_deviation_threshold)
	: record(the_record), quantity(the_quantity), deviation_threshold (the_deviation_threshold), calculation_time(false), total_time(false), all_scan(), countdown(NOT_SET), max_seconds(NOT_SET), p_base(0), p_state_start(0)
	{	};



/** initiate recursion **/
bool Engine::scan (
	vector<BaseData>& bases,
	const REAL contribution_threshold,
	const REAL threshold_factor
	)
{
	const char* here = "scan";
	if (threshold_factor >= 1.0) throw Exception (here, exc_ThresholdOne);
	if (contribution_threshold >= 1.0) throw Exception (here, exc_ThresholdOne);

	// overall stopwatch is running.
	total_time.start();

	// read interval file, if any.
	record.prog("reading interval file");
	const string file_name = quantity.fileName("int");
	all_scan.read(file_name);

	const Chain* p_chain = quantity.pLeftState()->p_chain;
/*
	// initial scan: scan only an approximately fixed (MAX_INITIAL) number of states on each base to determine the relative order
	for (int i=0; i< bases.size(); ++i) {
		// check constraints
		if (max_seconds != NOT_SET && calculation_time.seconds() > max_seconds) break;


		// skip initial scan if done something before on this base
		if (all_scan.numberAccepted(bases[i])) continue;
		record.prog()<< "initial scan of base "<< name(bases[i]) <<endl;

		// set a max number of states to scan
		countdown = MAX_INITIAL;
		p_base = newBase(p_chain, bases[i]);

		/// SKIP BASES WITH MORE THAN 2 TWO_STRINGS (XXX SPECIFIC!)
		if ((p_base->numberStringsOfType(1) > 1) || (p_base->numberTypes()>2)) {
			record.prog()<<"skipping base "<<name(bases[i])<<endl;
			delete p_base;
			continue;
		}
		///

		// do the initial scan
		p_state_start = 0;
		scanRecursive(1.0);
		delete p_state_start;
		delete p_base;

		// merge the results into the overall intervals list, and clear for next round
		all_scan.merge (last_scan);
		last_scan.clear();

 		record.prog()<< "base "<< name(bases[i]) <<SEP<< all_scan.summary(bases[i]) <<endl;
		record.prog()<< "cumulative "<< all_scan.summary() <<SEP<< calculation_time.humanReadable() <<endl;
	}
	record.prog("end of initial scans");
	record.prog()<<"cumulative "<< all_scan.summary() <<endl;
*/
	// no maximum number of states to scan
// 	countdown=NOT_SET;

	// sort bases by number of holes
	sortBases(bases, 0, bases.size()-1, BaseComparator(&all_scan));

	// steadily decrease threshold in steps of threshold_factor (<1.0)
	for (REAL threshold = threshold_factor; threshold >= contribution_threshold; threshold*=threshold_factor) {
		record.prog()<<"lowering threshold to "<< threshold <<endl;
		for (int i=0; i< bases.size(); ++i) {
			// check constraints
			if (max_seconds != NOT_SET && total_time.seconds() > max_seconds) {
				record.prog("total time limit reached.");
				break;
			}

// 			// have we checked on this base before and does it not yield anything significant? then skip.
			long long int done_before = all_scan.numberAccepted(bases[i]);
// 			if ( done_before && (!all_scan.numberCalculated(bases[i]) || all_scan.averageContribution(bases[i]) < threshold) ) continue;
			p_base = newBase(p_chain, bases[i]);

			/// SKIP BASES WITH MORE THAN 2 TWO_STRINGS (XXX SPECIFIC!)
/*			if ((p_base->numberStringsOfType(1) > 1) || (p_base->numberTypes()>2)) {
				record.prog()<<"skipping base "<<name(bases[i])<<endl;
				delete p_base;
				continue;
			}
*/
			///

			// if the base has been exhausted, don't bother.
			if (done_before == p_base->limId().back()) {
				record.prog()<<"all is done for base "<<name(bases[i])<<endl;
				delete p_base;
				continue;
			}
			record.prog()<<"scanning base "<< name(bases[i]) <<endl;

			// initiate recursive scan
			p_state_start = 0;
			scanRecursive(threshold);
			delete p_state_start;
			delete p_base;

			record.prog()<< "base "<<name(bases[i]) <<SEP<< last_scan.summary() <<endl;

			// did our most recent scan produce unsatisfactory results?
			bool last_under_threshold = last_scan.averageContribution() < threshold;

			// collect our results
			all_scan.merge (last_scan);
			last_scan.clear();

			record.prog()<< "cumulative "<< all_scan.summary() <<endl;
			record.prog()<< "calculation time "<< calculation_time.humanReadable() <<endl;
			record.prog()<< "total time "<< total_time.humanReadable() <<endl;

			// is this the first tyime we check this base and is it unsatisfying? skip all the rest.
			// NOTE that we assume that those done before are ordered before those not done before.
			if (!done_before && last_under_threshold) break;
		}
		record.prog()<< "threshold "<<threshold <<": cumulative "<< all_scan.summary() <<endl;
		record.prog()<< "calculation time "<< calculation_time.humanReadable() <<endl;
		record.prog()<< "total time "<< total_time.humanReadable() <<endl;

	}
	record.prog("writing interval file");
	all_scan.write(file_name);

	total_time.stop();
}



/** recursive part of scanning algorithm **/
void Engine::scanRecursive  (
	REAL contribution_threshold,
	int current_sector, int current_particle
	)
{
	const char* here = "scanRecursive";

	// nothing given: start at highest possible
	if (current_sector == NOT_SET) current_sector = p_base->numberSectors()-1; // number sectors must be higher than one... otherwise this base is bollox
	if (current_particle == NOT_SET) current_particle = max(0, -1 + (current_sector?p_base->number_particles_sector[current_sector]:p_base->numberHoles()) ) ;
	// if no particles in current sector, current_particle is zero. We should check against dereferencing this index below.

	/* negative indices: */
	//we shouldn't arrive here. something's wrong. bork.
	if ( (current_sector<0) || (current_particle<0) )
		throw Exception (here, exc_NegativeIndex);
	if ( (current_sector >= p_base->numberSectors())
		|| (current_particle && current_particle >= (current_sector?p_base->number_particles_sector[current_sector]:p_base->numberHoles()) ) )
		throw Exception (here, exc_IndexTooHigh);

	// ** the logic **
	// we have an current_sector (current_sector), pointing at which tableau we currently operate
	// and an current_particle (current_particle), pointing at which row we operate.
	// we check if the average contribution of whatever we've calculated 'under' this row is worthwhile;
	// * if so, we loop over the row and recursively call scanRecursive with a lower current_particle,
	//     or a lower current_sector and max current_particle if current_particle happens to be zero already.
	//     If everything remains worthwhile, we end up at indices zero and a scan will be done.
	// * if not so, we quit. this means we're going back to a point higher up in the recursion

	// get the shifts for the minimum and limit state we're considering
	// to find the maximum: in current shift, set rows above (particle indices below) the current row to maximum:

	// to find the minimum: in current shift, set rows above the current row to zero:
	vector<Young> start_shifts;
	if (!p_state_start) {
		/// TODO: this is initialisation code... move out? esp. creation of p_state?

		start_shifts = p_base->shifts(0);
		if (start_shifts[current_sector].height) {
			for (int row = 0; row <= current_particle; ++row)
				start_shifts[current_sector].writeAt(row) =
					(current_particle==start_shifts[current_sector].height-1)?0:start_shifts[current_sector].readAt(current_particle+1);
		}
		// set shifts below the current shift to zero
		for (int young=0; young < current_sector; ++young) start_shifts[young].setEmpty();
		p_state_start = newState(p_base, start_shifts);
	}
	else {
		start_shifts = p_state_start->shifts();
		if (start_shifts[current_sector].height) {
			for (int row = 0; row <= current_particle; ++row)
				start_shifts[current_sector].writeAt(row) =
					(current_particle==start_shifts[current_sector].height-1)?0:start_shifts[current_sector].readAt(current_particle+1);
		}
		p_state_start->setShift (current_sector, start_shifts[current_sector]);
	}

	vector<Young> stop_shifts = start_shifts;
	if (stop_shifts[current_sector].height)
		for (int row = 0; row <= current_particle; ++row)
			stop_shifts[current_sector].writeAt(row) = stop_shifts[current_sector].width;
	// set shifts below the current shift to full
	// setFull could be problematic for large numbers of holes?
	for (int young=0; young < current_sector; ++young)
		stop_shifts[young].setFull();
	long long int id_stop = p_base->id(stop_shifts) + 1;
	long long int id_start = p_state_start->id();
	Young shift = start_shifts[current_sector];

	/* all indices zero: */
	// we have arrived at the lowest level and must unfortunately start doing real work
	// we don't have to check averages here, since we haven't yet calculated these (or at least shouldn't have ;)
	// since these are the smallest building blocks
	if ( (0==current_sector) && (0==current_particle) ) {
		// innermost hole interval is atomic; hence, if one of its ids is included, it is fully included.
		if (!all_scan.includes (*p_base, id_start)) {

			scanInnermostHole (p_state_start, start_shifts[0].readAt(0), stop_shifts[0].readAt(0));

		}
		// subtract number of calculated states from countdown value. 0 signals a stop, negative is dangerous (NOT_SET is a negative number)
		if (countdown != NOT_SET) countdown = max(countdown-(stop_shifts[0].readAt(0)-start_shifts[0].readAt(0)+1), 0);
		return;
	}

	/* otherwise: */
	// we haven't yet arrived on the lowest level.
	// check if results from the past are encouraging
	// check if average contribution is high enough: if scans have been done and their contribution is too low, return.
	REAL average_contrib = all_scan.averageContribution(*p_base, id_start, id_stop);
	if ( finite(average_contrib) && (average_contrib < contribution_threshold) ) {
	 	return;
	}

	// find the right recursive call to maximise contributions.
	int new_current_sector = current_sector, new_current_particle = current_particle;
	if (0==current_particle) {
		// we can't lower the row; lower the tableau and set row to its max
		// (current_sector can't be zero, for then we would have been caught above)
		--new_current_sector;
		new_current_particle=NOT_SET; // implies max possible
	}
	else --new_current_particle;

	// NOTE that lower shifts have been set to zero above
	// find the highest id scanned before within the interval we're about to loop over.
	// (to check a bit further on if we've improved far enough)
	long long int highest_id_done = all_scan.highestId(*p_base, id_start, id_stop);
 	int highest_shift = (NO_ID==highest_id_done)? -1 : p_base->shifts(highest_id_done)[current_sector].readAt(current_particle);

	int blocks_calculated = 0;
	Interval last_diff = ZERO_INTERVAL;

	// loop over current row
	shift = start_shifts[current_sector];
	for (;;) {

		// set lower rows to their minimum length
		// (i.e. shift particles that further away from the middle by at least the shift of the current particle)
		// NOTE: a lower row is what you would draw higher in a Young diagram (yes, I know.)
		for (int lower_particle = 0; lower_particle < current_particle; ++lower_particle)
			shift.writeAt(lower_particle) = shift.readAt(current_particle);

		// change our running state
		p_state_start->setShift(current_sector, shift);
		// save summaries to be able to check contribution every other block
		Interval diff = last_scan.summary(*p_base);

		// recurse deeper.
		scanRecursive  (contribution_threshold, new_current_sector, new_current_particle);

		// check constraints (trackback recursion if so)
		if (max_seconds != NOT_SET && total_time.seconds() > max_seconds) break;
		/// TODO: why exactly zero?
		if (countdown==0) break;

		// check whether the contribution we just got was disappointing. we may want to break the loop and return to higher level.
		diff = last_scan.summary(*p_base) - diff;

		if (diff.number_calculated > 0) ++blocks_calculated;
		if (  ((blocks_calculated > 1)
			&& (diff.contribution+ last_diff.contribution)/(1.0*(diff.number_calculated+last_diff.number_calculated)) < contribution_threshold )

		   	|| (!diff.number_calculated && diff.number_accepted ) ) {
			// break only if we've calculated at least two blocks since the last run
			// (so that slots on both sides of the Fermi interval have been checked)
			if (shift.readAt(current_particle) - highest_shift >= 2 ) break;
		}

		last_diff = diff;

		if (shift.readAt(current_particle) >= stop_shifts[current_sector].readAt(current_particle)) break;
		++shift.writeAt(current_particle);
	}
}



/** innermost hole (terminating the recursion) **/
void Engine::scanInnermostHole (
	State* const p_state,
	const int start_hole,
	const int stop_hole)
{
	p_state->setFreeRapidities();
	Young shift = p_state->shifts()[0];
	long long int id_start=0;
	Interval done_interval = ZERO_INTERVAL;
	if (!shift.height) {
		p_state->setShift(0, shift);
		scanState(p_state, done_interval);
		done_interval.start = p_state->id();
	}
	else {
		// first scan the even shifts (so that we only move one quantum number a distance one at a time,
		// which, if we don't set free rapidities in the meantime, yields faster convergence.
		for (int hole_shift = start_hole; hole_shift <= stop_hole ; hole_shift+=2) {
			shift.writeAt(0) = hole_shift;
			p_state->setShift(0, shift);
			if (hole_shift==start_hole) done_interval.start = p_state->id();
			scanState(p_state, done_interval);
		}
		// then the odd shifts, the other way around
		for (int hole_shift = stop_hole - 1 + (start_hole+stop_hole)%2; hole_shift > start_hole; hole_shift-=2) {
			shift.writeAt(0) = hole_shift;
			p_state->setShift(0, shift);
			scanState(p_state, done_interval);
		}
	}
	done_interval.contribution /= quantity.maxSumHighestWeight ();
	done_interval.stop = done_interval.start+stop_hole-start_hole+1;
	last_scan.insert (*(p_state->p_base), done_interval);
}




void Engine::scanState(State* const p_state_in, Interval& scan_result)
{
	// weights for half Brillouin zone
	// add to statistics only if the contribution is in the correct half-Brillouin sector
	int weight = 0;
	if (p_state_in->mode()==0 || p_state_in->mode()==p_state_in->p_chain->length()/2) weight = 1;
	else if (p_state_in->mode()> 0  && p_state_in->mode() < p_state_in->p_chain->length()/2) weight = 2;
	else return;

	scan_result.number_accepted += weight;
	if (!p_state_in->admissible()) {
		// inadmissible state. log and drop.
		record.log (p_state_in) << ": inadmissible, q.n. " << p_state_in->quantum_number<<endl;
		return;
	}

	State* p_state = 0;
	calculation_time.start();
	try {
		// NOTE: WARNING! solve() may change the pointer itself to a newly allocated one !
		// this hack is necessary until we separate state id info from solving and roots.
		quantity.setRightState(p_state_in);
		p_state = solve(p_state_in, deviation_threshold, quantity.deviable());
		quantity.setRightState(p_state);

		// calculate form factor, cache it and write it away in a way we choose
		record(quantity);

		// update summaries
		scan_result.number_calculated += weight;
		if (finite(quantity.formFactor())) scan_result.contribution += quantity.formFactor() * REAL(weight);
		else record.log (p_state_in, "non-finite form factor");

	}
	catch (Exception exc) {
		// exception: log and move on.
		record.log (p_state_in, exc);
	};
	calculation_time.stop();
	// if we made something new, clean up.
	// deleting null pointers is harmless, double delete is not.
	if (p_state != p_state_in) delete p_state;
	quantity.setRightState(p_state_in);
}



