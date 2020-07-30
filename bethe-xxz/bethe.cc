
#include "bethe.h"

const char* exc_ZeroDenominator = "zero denominator"; //polint
const char* exc_NonFiniteArg = "non-finite argument"; // findContinuedFraction

const complex<REAL> I (0.0, 1.0);			// square root of -1
const REAL SQRT_PI=sqrt(PI);				
const REAL INFINITE = 1.0/0.0;				// I like compiler warnings
const REAL NOT_A_NUMBER = 0.0/0.0;

// global variable: null stream
nullstream CNULL;

/** continued fraction alorithm **/
// should perhaps be more properly termed findPartialQuotients
vector<int> findContinuedFraction(REAL REAL_number, REAL precision, int max_length)
{
	vector<int> continued_fraction;
	const char* here = "findContinuedFraction";
	
	// dumb cases
	if (REAL_number == 0.0) {
		continued_fraction.push_back(0);
		return continued_fraction;
	}
	if (!finite(REAL_number)) throw Exception (here, exc_NonFiniteArg);
	
	REAL remainder = 1.0/REAL_number; // the running variable
	int length = 0;
	REAL convergent = 0;

	while ((fabs(REAL_number - convergent) > precision) && (length < max_length)) {
		// calculate nth term in the continued fraction 
		continued_fraction.push_back (int(floor(remainder)));
		++length;
		remainder = 1.0/(remainder - floor(remainder));
		// nth convergent == continued fraction truncated at n
		// in other words, current approximation of REALNumber
		convergent = 0;
		for (int i = length-1; i >= 0; --i) convergent = 1.0/(continued_fraction[i] + convergent);
	}

	// should we throw an exception for this?
	// or just not give both precision and cutoff...
	// throw no exception. cutoff is more important than precision.
	//if (fabs(REAL_number - convergent) > precision) throw "findContinuedFraction: precision not reached due to cutoff";
	
	// if the last element is one, and there's more than one argument, remove it and add one to the second last.
	if ((continued_fraction.back() == 1) && (length>1)) {
		continued_fraction.pop_back(); 
		++continued_fraction.back(); 
	}
	
	return continued_fraction;
}


/** interpolate **/
// source: Numerical Recipes
 // Polynomial interpolation/extrapolation, NR page 113.
void polint (vector<REAL>& xa, vector<REAL>& ya, const REAL x, REAL& y, REAL& dy)
{
  int i, m, ns=0;
  REAL den, dif, dift, ho, hp, w;

  int n = xa.size();
  vector<REAL> c(n), d(n);
  dif = fabs(x-xa[0]);
  for (i = 0; i < n; i++) {
    if ((dift = fabs(x-xa[i])) < dif) {
      ns = i;
      dif = dift;
    }
    c[i] = ya[i];
    d[i] = ya[i];
  }
  y = ya[ns--];
  for (m = 1; m < n; m++) {
    for (i = 0; i < n-m; i++) {
      ho = xa[i] - x;
      hp = xa[i+m] - x;
      w = c[i+1] - d[i];
      if ((den = ho-hp) == 0.0) throw Exception ("polint", exc_ZeroDenominator);
      den = w/den;
      d[i] = hp * den;
      c[i] = ho * den;
    }
    y += (dy = (2 * (ns + 1) < (n-m) ? c[ns + 1] : d[ns--]));
  }
}


/** UNIX interface, generics **/

// FIXME: this doesn't work correctly if there are spaces in the file name...
vector<string> ls (const string description) 
{
	// do an ls command and write the output into a stringstream
	// okay, this is ugly. TODO: implement the right way using dirent.h
	stringstream response;	
	int c;
	FILE* f = popen(("ls " + description).c_str(), "r");
	while ((c = fgetc(f)) != EOF) response << (char)c;
	pclose(f);
	
	// get separate entries out of the stringstream into a vector.
	vector<string> listing;
	string entry;
	while (response >> entry) listing.push_back(entry);
	return listing;	
}

	
bool exists (const string file_name)
{
	FILE* result = fopen(file_name.c_str(), "r");
	if (!result) return false;
	
	fclose(result);
	return true;
}

string uniqueName (const string try_name)
{
	string the_name = try_name; 
	int i = 0; 
	while (exists(the_name)) {
		stringstream number;
		number << i++;
		the_name = try_name + '.' + number.str();
	}
	return the_name;
}

// tokenise a string (PHP explode)
// modified from http://oopweb.com/CPP/Documents/CPPHOWTO/Volume/C++Programming-HOWTO-7.html
vector<string> explode (const string& str, const string& delimiters)
{
	vector<string> tokens (0);
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
    return tokens;
}

/** class Stopwatch **/
/* just to get rid of the ugly C UNIX interface  */

string humanReadable(double secs)
{
	double days = floor ( secs / 86400.0 );
	secs -= days * 86400.0;
	double hours = floor( secs / 3600.0 ); 
	secs -= hours * 3600.0;
	double minutes = floor ( secs / 60.0 );
	secs -= minutes * 60.0;

	stringstream result;
	result << int(days) << "d " << int(hours) << "h " << int(minutes) <<"min ";
	result.precision(4);
	result << secs << "s ";
	return result.str();
}
 

Stopwatch::Stopwatch (bool running) : run(running) 
{
	p_usage = new rusage();
	reset();
	stop_u = start_u;
	stop_s = start_s;
}
	
Stopwatch::~Stopwatch(void) 
{
	delete p_usage;
}

void Stopwatch::start(void) 
{
	getrusage(RUSAGE_SELF, p_usage);
	// make up for time wasted
	REAL wasted_u = -stop_u + p_usage->ru_utime.tv_sec + (p_usage->ru_utime.tv_usec)/1.0e+6;
	REAL wasted_s = -stop_s + p_usage->ru_stime.tv_sec + (p_usage->ru_stime.tv_usec)/1.0e+6;
	start_u += wasted_u;
	start_s += wasted_s;
	lap_u += wasted_u;
	lap_s += wasted_s;
	// up and running
	run = true; 
}

void Stopwatch::stop(void) 
{
	getrusage(RUSAGE_SELF, p_usage);
	stop_u = p_usage->ru_utime.tv_sec + (p_usage->ru_utime.tv_usec)/1.0e+6;
	stop_s = p_usage->ru_stime.tv_sec + (p_usage->ru_stime.tv_usec)/1.0e+6;
	run = false;
}


void Stopwatch::reset(void) {
	getrusage(RUSAGE_SELF, p_usage);
	start_u = p_usage->ru_utime.tv_sec + (p_usage->ru_utime.tv_usec)/1.0e+6;
	start_s = p_usage->ru_stime.tv_sec + (p_usage->ru_stime.tv_usec)/1.0e+6;
	lap_u = start_u;
	lap_s = start_s;
}

double Stopwatch::getUse(void)  const
{  
	getrusage(RUSAGE_SELF, p_usage);
	if (run) return  - start_u + p_usage->ru_utime.tv_sec + (p_usage->ru_utime.tv_usec)/1.0e+6;	
	else return  stop_u - start_u;
}

double Stopwatch::getSys(void)  const
{  
	getrusage(RUSAGE_SELF, p_usage);
	if (run) return  - start_s + p_usage->ru_stime.tv_sec + (p_usage->ru_stime.tv_usec)/1.0e+6;	 
	else return stop_s - start_s;
}

double Stopwatch::lapUse(void)  
{  
	getrusage(RUSAGE_SELF, p_usage);
	double last_lap_u = lap_u;
	lap_u = p_usage->ru_utime.tv_sec + (p_usage->ru_utime.tv_usec)/1.0e+6;
	if (!run) lap_u = stop_u;
	return lap_u - last_lap_u;	 
}

double Stopwatch::lapSys(void)  
{  
	getrusage(RUSAGE_SELF, p_usage);
	double last_lap_s = lap_s;
	lap_s = p_usage->ru_stime.tv_sec + (p_usage->ru_stime.tv_usec)/1.0e+6;
	if (!run) lap_s = stop_s;
	return lap_s - last_lap_s;	 
}

