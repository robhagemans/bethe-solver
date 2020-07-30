

/** recursive scan over fixed base **/
		
		
void scanRecursiveTEST  (
	Engine& engine,
	const ScanIntervals& all_scan_interval,	// all the ones done before
	const ScanIntervals& run_scan_interval,	
	ScanIntervals* p_last_scan_interval,	// expects the ones done the last run,
	REAL contribution_threshold,
	bool extend,
	State* p_state_start,
	int current_sector,
	int current_particle)
{
	const char* here = "scanRecursive";
		

	/* checks on input */
	
	if (!p_state_start) throw Exception (here, exc_NullPointer);
	// for convenience
	Base* p_base = &p_state_start->base;
	
	// nothing given: start at highest possible
	// number sectors must be higher than one... otherwise this base is bollox
	if (current_sector == NOT_SET) current_sector = p_base->numberSectors()-1; 
	// if no particles in current sector, current_particle is zero. We should check against dereferencing this index below.
	if (current_particle == NOT_SET) current_particle = max(0, -1 + (current_sector?p_base->number_particles_sector[current_sector]:p_base->numberHoles()) ) ;

	// negative indices are illegal. bork.
	if ( (current_sector<0) || (current_particle<0) ) 
		throw Exception (here, exc_NegativeIndex);
		
	// check if indices aren't too high for this base
	if ( (current_sector >= p_base->numberSectors()) 
		|| (current_particle && current_particle >= (current_sector?p_base->number_particles_sector[current_sector]:p_base->numberHoles()) ) ) 
		throw Exception (here, exc_IndexTooHigh);
	
	
	// ** the logic **
	// we have an current_sector (current_sector), pointing at which tableau we currently operate
	// and an current_particle (current_particle), pointing at which row we operate.
	// we check if the average contribution of whatever we've calculated 'under' this row is worthwhile; 
	// * if so, we loop over the row and recursively call scanRecursive with a lower current_particle, 
	//      or a lower current_sector and max current_particle if current_particle happens to be zero already. 
	//      If everything remains worthwhile, we end up at indices zero and a scan will be done.
	// * if not so, we quit. this means we're going back to a point higher up in the recursion
	
	
	// get the shifts for the minimum and limit state we're considering
	// to find the maximum: in current shift, set rows above (particle indices below) the current row to maximum:
	
	// to find the minimum: in current shift, set rows above the current row to zero:
	// set shifts below the current shift to zero
	vector<Young> start_shifts = p_state_start->shifts();
	if (start_shifts[current_sector].height) 
		for (int row = 0; row <= current_particle; ++row) 
			start_shifts[current_sector].writeAt(row) = 
				(current_particle==start_shifts[current_sector].height-1)?0:start_shifts[current_sector].readAt(current_particle+1);
	for (int young=0; young < current_sector; ++young) start_shifts[young].setEmpty();
	p_state_start->setShift (current_sector, start_shifts[current_sector]);
	
	// find the sifts on the
	// set shifts below the current shift to full --- setFull could be problematic for large numbers of holes?
	vector<Young> stop_shifts = start_shifts; 		
	if (stop_shifts[current_sector].height) 
		for (int row = 0; row <= current_particle; ++row) 
			stop_shifts[current_sector].writeAt(row) = stop_shifts[current_sector].width;

	for (int young=0; young < current_sector; ++young) 
		stop_shifts[young].setFull();
	
	long long int id_stop = p_base->id(stop_shifts) + 1;
	long long int id_start = p_state_start->id();

	/* all indices zero: */
	// we have arrived at the lowest level and must unfortunately start doing real work
	// we don't have to check averages here, since we haven't yet calculated these (or at least shouldn't have ;)
	// since these are the smallest building blocks
	if ( (0==current_sector) && (0==current_particle) ) {
		// innermost hole interval is atomic; hence, if one of its ids is included, it is fully included.
		if ( !all_scan_interval.includes (p_base, id_start) &&
			( !run_scan_interval.includes(p_base, id_start)) ) {
			// innermost scan function: may be remote;
			// may return before calculations are done.  results collected elsewhere.
			engine.startScan (p_state_start, start_shifts[0].readAt(0), stop_shifts[0].readAt(0));
		}
		return;
	}

	/* otherwise: */
	// we haven't yet arrived on the lowest level. 
	// check if results from the past are encouraging
	// check if average contribution is high enough: if scans have been done and their contribution is too low, return.
	REAL average_contrib = all_scan_interval.averageContribution(p_base, id_start, id_stop);
	if ( finite(average_contrib) && (average_contrib < contribution_threshold) ) return;

	// find the next sector and particle
	int new_current_sector = current_sector, new_current_particle = current_particle;
	if (0==current_particle) {
		// we can't lower the row; lower the tableau and set row to its max
		// (current_sector can't be zero, for then we would have been caught above)
		--new_current_sector;
		new_current_particle=NOT_SET; // implies max possible
	}
	else --new_current_particle;	

	// NOTE that lower shifts in  start_shifts / p_state_start.shifts()  have been set to zero above 

	
	// two types of scan; extend scan goes two boxes beyond last; sniff scan checks last to see if worthwhile
	int last_box;
	if (extend) {
		
		/** extend scan **/
		// find the highest id scanned before within the interval we're about to loop over.
		long long int highest_id_done = all_scan_interval.highestId(p_base, id_start, id_stop);
		if (NO_ID!=highest_id_done) last_box = p_base->shifts(highest_id_done)[current_sector].readAt(current_particle) +2;
		else last_box = 1; // nothing done: do 0 and 1
	}
	else {
		
		/** sniff scan **/
		// sniff scan: find the furthest we went earlier *this threshold/base run* 
	 	long long int highest_id_done = run_scan_interval.highestId(p_base, id_start, id_stop);
	 	if (NO_ID!=highest_id_done) last_box = p_base->shifts(highest_id_done)[current_sector].readAt(current_particle);
		else last_box = -1; 
		
		// if we have already calculated the last box on this row before, there's nothing to do anymore.
		if (last_box >= stop_shifts[current_sector].readAt(current_particle)) return;
		
		if (last_box == -1) {
			/// this almost works .	
			/// TODO: figure out why not exactly	
			
			// last_box = -1: this row hasn't been scanned before, 
			// so we can't check the last calculated box on the row
			// instead, do an extension of two boxes 0 and 1. (as in the initial extend scan)
			extend=true;
			last_box=1;
		}
		// if we're doing a sniff scan, check if the last scanned before was worthwhile before doing the next
		else {
			// NOTE that it is essential that last_box+1 is still in the allowed interval (checked above)
			// notably, last_box may not be <0 or >= stop_shifts[current_sector][current_particle]
			
			vector<Young> shift_temp = start_shifts;
			for (int lower_particle = 0; lower_particle <= current_particle; ++lower_particle) 
				shift_temp[current_sector].writeAt(lower_particle) = last_box; 
			long long int start_id_last_box = p_base->id(shift_temp);
			
			for (int lower_particle = 0; lower_particle <= current_particle; ++lower_particle) 
				shift_temp[current_sector].writeAt(lower_particle) = last_box+1; 
			long long int stop_id_last_box = p_base->id(shift_temp) +1;	
	
			if (stop_id_last_box==start_id_last_box) throw Exception(here, exc_LastBoxEmpty);
	
			// should we go one box further?
			Interval last_box_summary = run_scan_interval.summary(p_base, start_id_last_box, stop_id_last_box);
			if (!( ((last_box_summary.number_calculated > 0) 
				&& (last_box_summary.contribution)/(1.0*last_box_summary.number_calculated) < contribution_threshold )
				|| (!last_box_summary.number_calculated && last_box_summary.number_accepted ) ) ) 
						++last_box;
		}
	}
		
	// don't go further than we can
	last_box = min(last_box, stop_shifts[current_sector].readAt(current_particle));
	last_box = max(last_box, start_shifts[current_sector].readAt(current_particle));
	
	// loop over current row
	for (	Young shift = start_shifts[current_sector];
			shift.readAt(current_particle) <= last_box; 
			++shift.writeAt(current_particle)) {

		// set lower rows to their minimum length
		// (i.e. shift particles that further away from the middle by at least the shift of the current particle)
		// NOTE: a lower row is what you would draw higher in a Young diagram (yes, I know.)
		for (int lower_particle = 0; lower_particle < current_particle; ++lower_particle) 
			shift.writeAt(lower_particle) = shift.readAt(current_particle); 
			
		// first copy our state to make a running state (scanRecursive will alter it, so we can't use our original)
		// then change to our running state
		State* p_state = copy (p_state_start);		
		p_state->setShift(current_sector, shift);
		
		// recurse deeper.
		scanRecursiveTEST  (engine, all_scan_interval, run_scan_interval, p_last_scan_interval, 
			 	contribution_threshold, extend, p_state, new_current_sector, new_current_particle);
			 	
		// delete our used copy
		delete p_state;
	}
}




void scanRecursiveTEST (		
	Engine& engine,
	const ScanIntervals& all_scan_interval,	// all the ones done before
	const ScanIntervals& run_scan_interval,	
	ScanIntervals* p_last_scan_interval,	// expects the ones done the last run,
	REAL contribution_threshold,
	bool extend,
	Base* p_base)
{
	// initial call (no start state given, set ground state)
	State* p_state_start = newState(*p_base, 0);
	scanRecursiveTEST (engine, all_scan_interval, run_scan_interval, p_last_scan_interval, contribution_threshold, extend, p_state_start, NOT_SET, NOT_SET);	
	delete p_state_start;
}

