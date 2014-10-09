#ifndef STRIP_H
#define STRIP_H

#include <vector>
#include <algorithm>
#include <string>
using namespace std;

#include "basics.h"
#include "exception.h"
#include "configuration.h"


extern const string exc_OutOfRange;
extern const string exc_BadAssign;


template<class number> 
class Strip {
public:

	// construct from Configuration
	Strip (const Configuration& new_base, const number initial_value=0);
	// copy constructor
	Strip (const Strip& original);	
	// destructor
	~Strip();	
	// assignment			
	Strip<number>& operator= (const Strip<number>& original); // assignment

    Strip<number>& set(const int j, const int alpha, const number value);
	
	inline number element(int i) const 
	{ 	return elements_[i]; }
	
	inline number& element(int i) 
	{ 	return elements_[i];	}
	
	inline int index(const int j, const int alpha) const 
	{ 	return alpha + index_table_[j]; };
	
	inline number& operator() (const int j, const int alpha) 
	{ 	return elements_[index(j, alpha)];	};
	
	inline number operator() (const int j, const int alpha) const 
	{  	return elements_[index(j, alpha)]; 	};
	
	// number of string types
	inline int numberTypes(void) const  
	{ 	return number_types_; 	};
	
	// number of strings of jth string type. should be equal to base.numberStringsOfType()
	inline int numberStringsOfType (const int j) const 	
	{ 	return (j<number_types_)?(index_table_[j+1] - index_table_[j]):0; 	};
	
	inline int numberElements(void) const  
	{ 	return (index_table_[number_types_] - index_table_[0]);	};

protected:
	number* elements_;
	int* index_table_;
	int number_types_;
	//const Configuration config_;
	
	template<class Compare>
	friend void sort(Strip<number>& strip_to_sort, Compare compare);
};




template<class number, class Compare> 
void sort(Strip<number>& strip_to_sort, Compare compare);


	


template<class number> 
class FullStrip {
public:
	// construct from Configuration
	FullStrip (const Configuration& new_base, const number initial_value=0);
	// copy constructor
	FullStrip (const FullStrip& original);	
	// destructor
	~FullStrip();
	// assignment				
	FullStrip<number>& operator= (const FullStrip<number>& original); 
	
	// index triplet to index 
	inline int index(const int j, const int alpha, const int a) const 
	{  	return  a + alpha*string_length_[j] + index_table_[j]; };
	
	// element extraction by index triplet (write)
	inline number& operator() (const int j, const int alpha, const int a) 
	{  	return elements_[index(j, alpha,a)]; };

	// element extraction by index triplet (read)
	inline number operator() (const int j, const int alpha, const int a) const 
	{  	return elements_[index(j, alpha,a)]; };
	
	// number of string types
	inline int numberTypes() const  
	{ 	return number_types_; };
	
	// number of strings of jth string type. should be equal to base.numberStringsOfType()
	inline int numberStringsOfType (const int j) const  
	{ 	return (j<number_types_)?(index_table_[j+1] - index_table_[j]):0;  };
	
	inline int numberElements() const  
	{ 	return (index_table_[number_types_] - index_table_[0]);};
	
protected:
	number* elements_;
	int* index_table_;
	int number_types_;
	vector<int> string_length_;
};






template<class number> 
Strip<number>::Strip (const Configuration& config, const number initial_value) 
	: number_types_(config.numberTypes())
{
	int number_elements = 0;
	for (int j=0; j< number_types_; ++j) {
		number_elements += config.numberStringsOfType(j);
	}

	elements_ = new number[number_elements];
	index_table_ = new int[config.numberTypes()+1];
	index_table_[0]=0;
	// look out! the last pointer points one past the allocated space. // useful to get size info. 
	for (int i = 1; i <= number_types_; ++i) 
		index_table_[i] = index_table_[i-1] + config.numberStringsOfType(i-1);
	// zero out.
	for (int i = 0; i < number_elements; ++i) 
		elements_[i] = initial_value;
}




// copy constructor
template<class number> 
Strip<number>::Strip (const Strip& original): number_types_(original.numberTypes())
{
	index_table_ = new int[number_types_+1];
	elements_ = new number[original.numberElements()];
	for (int i = 0; i <= number_types_; ++i) 
		index_table_[i] = original.index_table_[i];
	for (int i = 0; i < original.numberElements(); ++i) 
		elements_[i] = original.elements_[i];
}



// assignment
template<class number> 
Strip<number>& Strip<number>::operator= (const Strip<number>& original)
{
	if (numberElements() != original.numberElements()) 
		throw Exception("Strip::operator=", exc_BadAssign, "(number elements)");
	if (number_types_ != original.numberTypes()) 
		throw Exception("Strip::operator=", exc_BadAssign, "(number types)");

	if (this != &original) {
		for (int i = 0; i <= number_types_; ++i) {
			if (index_table_[i] != original.index_table_[i]) 
				throw Exception("Strip::operator=", exc_BadAssign, "(index table)");
		}
		// copy elements
		for (int i = 0; i < original.numberElements(); ++i) 
			elements_[i] = original.elements_[i];
	}

}


// set a value in the strip - chainable method
template<class number>
Strip<number>& Strip<number>::set(const int j, const int alpha, const number value) 
{
    elements_[index(j, alpha)] = value;
    return *this;
}

// destructor
template<class number>
Strip<number>::~Strip() { 
	delete[] elements_;
	delete[] index_table_;
}




// outputting a Strip
template<class number>
ostream& operator<< (ostream& stream, const Strip<number>& out_to_be_put) {
	
	if (out_to_be_put.numberTypes() > 0) 
		for (int i=0; i< out_to_be_put.numberTypes(); ++i) {
		    
			for (int j=0; j < out_to_be_put.numberStringsOfType(i); ++j) {
 				stream << ' ' << out_to_be_put(i, j);
			}
			if (out_to_be_put.numberStringsOfType(i)>0) stream << ' ';
			stream << ";";
		}
	return stream;
}

// sorting a strip
template<class number, class Compare>
void sort(Strip<number>& strip_to_sort, Compare compare)
{
	for (int i=0; i<strip_to_sort.numberTypes(); ++i) 
		sort (strip_to_sort.elements_ + strip_to_sort.index_table_[i], strip_to_sort.elements_+strip_to_sort.index_table_[i+1], compare);
}





/* 
/* FullStrip
 */
 

template<class number> 
FullStrip<number>::FullStrip (const Configuration& config, const number initial_value) 
	:  number_types_(config.numberTypes()), string_length_(config.stringLengths())
{

	int number_elements = config.numberRoots();
	
	elements_ = new number[number_elements];
	index_table_ = new int[number_types_+1];
	
	index_table_[0]=0;
	// look out! the last pointer points one past the allocated space. // useful to get size info. 
	for (int i = 1; i <= number_types_; ++i) 
		index_table_[i] = index_table_[i-1] + config.numberStringsOfType(i-1) * string_length_[i-1];
	
	// zero out.
	for (int i = 0; i < number_elements; ++i) 
		elements_[i] = initial_value;
}


// copy constructor
template<class number> 
FullStrip<number>::FullStrip (const FullStrip& original)
	: number_types_(original.numberTypes()), string_length_(original.string_length_)
{
	elements_ = new number[original.numberElements()];
	index_table_ = new int[number_types_+1];
	
	// copy elements. 
	for (int i = 0; i < original.numberElements(); ++i) 
		elements_[i] = original.elements_[i];
	// copy index table
	for (int i = 0; i <= number_types_; ++i) 
		index_table_[i] =  original.index_table_[i];
}

// assignment
template<class number> 
FullStrip<number>& FullStrip<number>::operator= (const FullStrip<number>& original)
{
	if (this != &original) {
		// check index table
		for (int i = 0; i <= number_types_; ++i) {
			if (index_table_[i] != original.index_table_[i]) 
				throw Exception("FullStrip::operator=", exc_BadAssign, " (index table)");
		}
		// copy elements
		for (int i = 0; i < original.numberElements(); ++i) {
			elements_[i] = original.elements_[i];
		}
	}
}


// destructor
template<class number>
FullStrip<number>::~FullStrip() 
{
	delete[] elements_;
	delete[] index_table_; 
}





#endif

