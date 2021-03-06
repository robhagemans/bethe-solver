#ifndef STRIP_H
#define STRIP_H

#include <vector>
#include <algorithm>
using namespace std;
#include "bethe.h"
#include "base.h"

// use extremely range-checked versions to get more useful errors than "memory corruption. abort."
// #define STALINIST_RANGE_CHECK

// this has to do with the fact that GCC has to include the template definition in each object file
// I should have a .cc file specially for these constants.
#define exc_OutOfRangeASDEFINEDINTHISFILE "index out of range"
#define exc_BadAssignASDEFINEDINTHISFILE "bad assignment"

#ifndef STALINIST_RANGE_CHECK

template<class number>
class Strip {
public:
	number* elements;
	int* index_table;
	int number_types;
	const Base* p_base;
public:
	Strip (const Base* const new_p_base, const number initial_value=0);
	Strip (const Strip& original);	// copy constructor
	~Strip();
	Strip<number>& operator= (const Strip<number>& original); // assignment

	inline number element(int i) const { return elements[i]; }
	inline number& element(int i) { return elements[i];	}
	inline int index(const int j, const int alpha) const { 	return alpha + index_table[j]; };
	inline number& operator() (const int j, const int alpha) { return elements[index(j, alpha)];};
	inline number operator() (const int j, const int alpha) const {  return elements[index(j, alpha)]; };

	// number of string types
	inline int numberTypes(void) const  { return number_types; };

	// number of strings of jth string type. should be equal to base.numberStringsOfType()
	inline int numberStringsOfType (const int j) const 	{ return (j<number_types)?(index_table[j+1] - index_table[j]):0; };

	inline int numberElements(void) const  { return (index_table[number_types] - index_table[0]);	};
};


template<class number, class Compare>
	void sort(Strip<number>& strip_to_sort, Compare compare);



template<class number>
class FullStrip {
public:
	number* elements;
	int* index_table;
	int number_types;
	const Base* p_base;
public:
	// construct from Base
	FullStrip (const Base* const new_p_base, const number initial_value=0);
	// copy constructor
	FullStrip (const FullStrip& original);
	// destructor
	~FullStrip();
	// assignment
	FullStrip<number>& operator= (const FullStrip<number>& original);

	// index triplet to index
	inline int index(const int j, const int alpha, const int a) const {  return  a + alpha*p_base->p_chain->stringLength(j) + index_table[j]; };

	// element extraction by index triplet (write)
	inline number& operator() (const int j, const int alpha, const int a) {  return elements[index(j, alpha,a)]; };

	// element extraction by index triplet (read)
	inline number operator() (const int j, const int alpha, const int a) const {  return elements[index(j, alpha,a)]; };

	// number of string types
	inline int numberTypes(void) const  { return number_types; };

	// number of strings of jth string type. should be equal to base.numberStringsOfType()
	inline int numberStringsOfType (const int j) const  { 	return (j<number_types)?(index_table[j+1] - index_table[j]):0;  };

	inline int numberElements(void) const  { return (index_table[number_types] - index_table[0]);};
};

#else

/** ********************** ************************* *************************** ************************** **/
//#ifdef STALINIST_RANGE_CHECK

/* Failsafe versions:
/* these are totalitarian-style range-chacked versions of the above for debugging purposes (death to memory corruptors!)
*/


template<class number>
class Strip {
public:
	vector<number> elements;
	vector<int> index_table;
	int number_types;
	Base* p_base;
public:
	Strip (Base* new_p_base, const number initial_value=0);
	Strip (const Strip& original);	// copy constructor
	~Strip();
	Strip<number>& operator= (const Strip<number>& original); // assignment

	inline number element(int i) const
	{
		if (i<0 || i>elements.size()) throw Exception("Strip::element() const", exc_OutOfRangeASDEFINEDINTHISFILE);
		return elements[i];
	}

	inline number& element(int i)
	{
		if (i<0 || i>elements.size()) throw Exception("Strip::element()", exc_OutOfRangeASDEFINEDINTHISFILE);
		return elements[i];
	}

	inline int index(const int j, const int alpha) const
	{
		const char* here = "Strip::index()";
		try {
			if (j>=numberTypes() || j<0) throw Exception(here, exc_OutOfRangeASDEFINEDINTHISFILE, "(j)");
			if (alpha>=p_base->numberStringsOfType(j) || alpha<0) throw Exception(here, exc_OutOfRangeASDEFINEDINTHISFILE, "(alpha)");

			int ze_index =  alpha + index_table.at(j);
			if (  (ze_index >= numberElements())  || ze_index < 0) throw Exception(here, exc_OutOfRangeASDEFINEDINTHISFILE, "(index)");
			return ze_index;
		}
		catch (std::out_of_range) {
			throw Exception(here, exc_OutOfRangeASDEFINEDINTHISFILE, "(index table)");
		}
	};

	inline number& operator() (const int j, const int alpha)
	{
		if (index(j,alpha) > numberElements()) throw Exception("Strip<number>::operator()", exc_OutOfRangeASDEFINEDINTHISFILE);
		try {
			return elements.at(index(j, alpha));
		}
		catch (std::out_of_range) {
			throw Exception("Strip<number>::operator()", exc_OutOfRangeASDEFINEDINTHISFILE);
		}
	};

	inline number operator() (const int j, const int alpha) const
	{
		if (index(j,alpha) > numberElements()) throw "range (read)";
	 	try {
	 		return elements.at(index(j, alpha));
	 	}
	 	catch (std::out_of_range) {
			throw Exception("Strip<number>::operator() const", exc_OutOfRangeASDEFINEDINTHISFILE);
		}
	};

	// number of string types
	inline int numberTypes(void) const
	{
		return number_types;
	};


	// number of strings of jth string type. should be equal to base.numberStringsOfType()
	inline int numberStringsOfType (const int j) const
	{
		try {
			return (j<number_types)?(index_table.at(j+1) - index_table.at(j)):0;
		}
		catch (std::out_of_range) {
			throw Exception("Strip<number>::numberStringsOfType()", exc_OutOfRangeASDEFINEDINTHISFILE, "(index table)");
		}
	};

	inline int numberElements(void) const
	{
		if (index_table.at(number_types) - index_table.at(0) != elements.size()) {
			throw Exception("Strip<number>::numberElements()", exc_OutOfRangeASDEFINEDINTHISFILE, "strip corrupted");
		}
		return (index_table.at(number_types) - index_table.at(0));

	 };
};


template<class number, class Compare>
	void sort(Strip<number>& strip_to_sort, Compare compare);



/* FullStrip, the neostalinist version
**/


template<class number>
class FullStrip {
public:
	vector<number> elements;
	vector<int> index_table;
	int number_types;
	Base* p_base;
public:
	// construct from Base
	FullStrip (Base* new_p_base, const number initial_value=0);
	// copy constructor
	FullStrip (const FullStrip& original);
	// destructor
	~FullStrip();
	// assignment
	FullStrip<number>& operator= (const FullStrip<number>& original);

	// index triplet to index
	inline int index(const int j, const int alpha, const int a) const
	{
		try {
			int ze_index =  a + alpha*p_base->chain.stringLength(j) + index_table.at(j);
			if (  (ze_index >= numberElements())  || ze_index < 0) throw Exception("FullStrip::index()", exc_OutOfRangeASDEFINEDINTHISFILE);
			return ze_index;
		}
		catch (std::out_of_range) {
			throw Exception("FullStrip::index()", exc_OutOfRangeASDEFINEDINTHISFILE, "(index table)");
		}
	};

	// element extraction by index triplet (write)
	inline number& operator() (const int j, const int alpha, const int a)
	{
		if (index(j,alpha,a) > numberElements()) throw Exception("FullStrip::operator()", exc_OutOfRangeASDEFINEDINTHISFILE);
		try {
			return elements.at(index(j, alpha,a));
		}
		catch (std::out_of_range) {
			throw Exception("FullStrip::operator()", exc_OutOfRangeASDEFINEDINTHISFILE, "(index table)");
		}
	};


	// element extraction by index triplet (read)
	inline number operator() (const int j, const int alpha, const int a) const
	{
		if (index(j,alpha,a) > numberElements()) throw Exception("FullStrip::operator() const", exc_OutOfRangeASDEFINEDINTHISFILE);
	 	try {
	 		return elements.at(index(j, alpha,a));
	 	}
	 	catch (std::out_of_range) {
			throw Exception("FullStrip<number>::operator() const", exc_OutOfRangeASDEFINEDINTHISFILE, "(index table)");
		}
	};


	// number of string types
	inline int numberTypes(void) const
	{
		return number_types;
	};


	// number of strings of jth string type. should be equal to base.numberStringsOfType()
	inline int numberStringsOfType (const int j) const
	{
		try {
			return (j<number_types)?(index_table.at(j+1) - index_table.at(j)):0;
		}
		catch (std::out_of_range) {
			throw Exception("FullStrip<number>::numberStringsOfType", exc_OutOfRangeASDEFINEDINTHISFILE, "(index table)");
		}

	};

	inline int numberElements(void) const
	{
		if (index_table.at(number_types) - index_table.at(0) != elements.size()) {
			throw Exception("FullStrip<number>::numberElements()", exc_OutOfRangeASDEFINEDINTHISFILE, "strip corrupted");
		}
		return (index_table.at(number_types) - index_table.at(0));
	};
};


#endif


// need to include template definitions in header...
#include "strip.cc"

#endif
