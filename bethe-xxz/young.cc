#include "young.h"

const char* exc_Overflow = "overflow"; // factorial(), choose()
const char* exc_InvalidID = "invalid id"; //young()
//const char* exc_IndexRange = "index out of range"; //young() // defined in matrix.h



/** factorial **/
long long int factorial (int n) 
{
	const char* here = "factorial";
	
	long long int fac = 1;
	for (int i=2; i <= n; ++i) {
		fac *= i;
		if (fac<0) throw Exception (here, exc_Overflow); 
	}
	return fac;
}

/** choose (binomial coefficient) **/
long long int choose (int upper, int lower)
{
	const char* here = "choose";
	if ((upper < 0)||(lower<0)) return 0;
	if (upper < lower) return 0;
	if (upper == lower) return 1;
	
	// don't do unnecessary work with overflow risk
	if (lower>(upper-lower)) lower = upper-lower;
	long long int num = 1;
	for (int i=1; i <= int(lower); ++i) {
		num *= upper- i +1;
		
		// keep numbers from overflowing while not losing exactness
		if (num%i) throw Exception (here, exc_Overflow);
			// they really should always be divisible!
		else num /=i;
		if (num<0) throw Exception (here, exc_Overflow);
	}
	return (num);
}




/** class Young **/

/// NOTE: I currently create them at their max_size. this is bad for large (e.g. N=8000) chains, the tableaux are just too big!

// construct young tableau (a vector of its row widths) from id number
Young::Young (const int max_width, const int max_height, const long long int id): rows(max_height) 
{ 
	setId(max_width, max_height, id); 
	if (!height) rows.resize(1,0);
} ;
	
// convert vector<int> to Young
Young::Young (const int max_width, const int max_height, const std::vector<int>& row_widths) : rows(max_height), width(max_width), height(max_height) 
{ 
	for (int i=0; i< row_widths.size(); ++i) rows[i] = row_widths[i];
	if (!height) rows.resize(1,0);
};

// copy
Young::Young (const Young& original) : rows(original.rows), width(original.width), height(original.height) 
{ 
	if (!height) rows.resize(1,0);
};
	

Young& Young::operator= (const Young& rhs)
{
	if (this != &rhs) {
		width = rhs.width;
		height = rhs.height;
		// assign the vector part. 
		rows = rhs.rows;
	}
	return *this;
}


int Young::readAt (const int index) const
{
	const char* here = "Young::read";
	// it is legal to ask for value(0) of an empty tableau. it equals zero.
	if (!height && !index) return 0;
	if ((index >= height) || (index < 0)) {
		throw Exception (here, exc_IndexRange);
	}
	// if we read beyond the implemented length, we get zero.
	// (however, we're now always max_height long...)
	if (index >= rows.size()) return 0;
 	return rows.at(index);
}


int& Young::writeAt (const int index) {
	const char* here = "int& Young::at";
	
	if ((index >= height) || (index < 0)) {	
		throw Exception (here, exc_IndexRange);
	}
	// we try to write to index bigger than we have yet, but legal, so resize
	if (index >= rows.size()) rows.resize(index+1, 0);
 	return rows.at(index);
}


void Young::setId (const int max_width, const int max_height, long long int id)
{
	const char* here = "Young::setId";
	width=max_width;
	height=max_height;
	rows.clear();
	

	if ((id<0) || (id>=choose(max_height+max_width, max_width))) {
		stringstream remark;
		remark << 0 <<" <= " << id <<" < " << choose(max_height+max_width, max_width);
		throw Exception (here, exc_InvalidID, remark.str());
	}
	if (!height) return;
	
	// for historic reasons I construct the transpose, then transpose
	int last_height = max_height;
	long long int last_chs=0;
	long long int chs = 1;
	for (int column=0; column < max_width; ++column) {
		int h = 0;
		for (h = 0; h <= last_height; ++h) {
			last_chs = chs;
			chs = choose (max_width -column + h, max_width-column);
			// checked access resizes vector if necessary //
			if (h>0) (writeAt(h-1))++; // tp[column] = h;
			if (id < chs) break;
		}
		id -= last_chs;
		last_height = h;
		if (last_height==0) break;
	}
}

long long int Young::id (void) const
{
	if (!height) return 0;
	if (!rows.size()) return 0;
	
 	// we'd like to know the heights
 	Young tp (height, width);
	for (int i=0; (i < rows.size()) && (i<height) && readAt(i); ++i) 
		for (int j=0; j < readAt(i); ++j) 
			++tp.writeAt(j); // Young::at() (and Young::operator[]) resize if necessary
	
	// why does this work? anyone?
	long long int id = 0;
	for (int column=0; column < width; ++column)
		id += choose (width - column + tp.readAt(column)-1, width - column);
		// this equals the number of tableaux of width (width-column) ==the width right of and including the current column
		// and height (tp[column]-1) == height of current column *minus one*.
		// in other words, at every column we add the number of diagrams that could possibly be 
		// fit in the space from at its column, to the right; and from above its top, to top.
		// e.g. at width=5, height=3,
		//   ***xx
		//   ***o
		//   **
		// we're considering the box at o. then we count everything possible that fits in xx (3 diagrams).
		// which is the number of diagrams that went before this one, at this stage: 
		// we count 0 ; 1 ; 2  ; 3 ;  4 ;  5
		//          . ; * ; ** ; * ; ** ; ** 
		//                       *   *    **
		
	return id;
}

Young& Young::transpose (void)
{
	// copy
	Young tp (*this);
	// transpose dimensions
	int was_height=height;
	height=width;
	width=was_height;
	// don't bother with empty tableaux
	if (!height) return (*this);
	if (!rows.size()) return (*this);
	// rebuild as transpose
	rows.clear();
	for (int i=0; (i < tp.rows.size()) && (i<tp.height) && tp.readAt(i) ; ++i) 
		for (int j=0; j < tp.readAt(i); ++j) 
			++writeAt(j);
	return (*this);
}

Young& Young::setFull (void)
{
	rows.clear();
	// fill (height) rows with the value (width)
	rows.resize(height, width);
	return (*this);
}

Young& Young::setEmpty (void)
{
	rows.clear();
	rows.resize(1,0);
	return (*this);
}


// outputting a tableau
ostream& operator<< (ostream& stream, Young out_to_be_put) {
	stream<<" w "<<out_to_be_put.width<<" h "<<out_to_be_put.height;
	stream << "[";
	if (out_to_be_put.height > 0) stream << out_to_be_put.readAt(0);
	if (out_to_be_put.height > 1) 
		for (int i=1; i< out_to_be_put.height; ++i)	
			stream << "; " << out_to_be_put.readAt(i);
	
	stream << "]";
	return stream;
}

/*
	visualise a tableau:
	
	for (int i=0; i<out_to_be_put.size(); ++i) {
		for (int j=0; j< out_to_be_put[i]; ++j)
			stream << "O";
		stream<<endl;
	}
*/

const vector<Young> NO_SHIFTS (0);
