#ifndef STRIP_DEF
#define STRIP_DEF

#include "strip.h"

/* class Strip: implementation
/*
/**/

#ifndef STALINIST_RANGE_CHECK

/*

//these lead to incomprehensible errors (currently, badalloc)
// hmm, that might be due to bad linking because of makefile error....


template<class number>
Strip<number>::Strip (Base* new_p_base, const number initial_value)
{
	dimension = new_p_base->numberTypes();

	contents = new number*[dimension+1];
	int number_elements = 0;
	for (int j = 0; j < dimension; ++j) number_elements += new_p_base->numberStringsOfType(j);

	contents[0] = new number[number_elements];

	for (int i = 1; i <= dimension; ++i) contents[i] = contents[i-1] + new_p_base->numberStringsOfType(i-1);
	// look out! the last pointer points one past the allocated space.
	// useful to get size info. maybe implement this differently.

	// zero out.
	for (int i = 0; i < number_elements; ++i) *(contents[0] + i) = initial_value;
}


// copy constructor
template<class number>
Strip<number>::Strip (const Strip& original): dimension(original.numberTypes())
{
	/// NOTE: watch this constructor! It's caused memory corruption before
	/// DO NOT use contents(new number* [dimension+1]) it the initialiser, it corrupts memory (or at least I think it's the cause...)

	contents= new number* [dimension+1];

	contents[0] = new number[original.numberElements()];
	// copy elements.
	for (int i = 0; i < original.numberElements(); ++i) *(contents[0] + i) = *(original.contents[0] + i);
	// 'copy' pointer table, but with pointers to the new array
	for (int i = 1; i <= dimension; ++i) contents[i] = contents[0] + (original.contents[i] - original.contents[0]);
}

// assignment
template<class number>
Strip<number>& Strip<number>::operator= (const Strip<number>& original)
{
	if (this != &original) {
		if (numberElements() != 	original.numberElements()) throw "bad assignment (#)";
		if (dimension != 	original.numberTypes()) throw "bad assignment (dimension)";
		// copy elements
		for (int i = 0; i < original.numberElements(); ++i) *(contents[0] + i) = *(original.contents[0] + i);
		// 'copy' pointer table, but with pointers to the new array
		/// don't do this for an assignment! the structures must be equal anyway...
// 		for (int i = 1; i <= dimension; ++i) contents[i] = contents[0] + (original[i] - original[0]);
	}
}



// destructor
template<class number>
Strip<number>::~Strip()
{
	if (numberElements()) delete[] contents[0];
	delete[] contents;
}


// template<class number>
// vector<int> Strip<number>::structure(void) const
// {
// 	vector<int> result_structure (dimension);
// 	for (int i = 0; i < dimension; ++i) result_structure[i] = (contents[i+1] - contents[i]);
// 	return result_structure;
// }



// outputting a Strip
template<class number>
ostream& operator<< (ostream& stream, const Strip<number>& out_to_be_put) {

	//stream << "(";
	if (out_to_be_put.numberTypes() > 0)
		for (int i=0; i< out_to_be_put.numberTypes(); ++i) {
			for (int j=0; j < out_to_be_put.numberElementsOfType(i); ++j)
				stream << out_to_be_put[i][j] << SEP;//", ";
			stream << ";" << SEP;//";  ";
		}
	//stream << ")";

	return stream;
}

// sorting a strip
template<class number, class Compare>
void sort(Strip<number>& strip_to_sort, Compare compare)
{
	for (int i=0; i<strip_to_sort.numberTypes(); ++i) sort (strip_to_sort.contents[i], strip_to_sort.contents[i+1], compare);
}

*/




template<class number>
Strip<number>::Strip (const Base* const new_p_base, const number initial_value)
// 	: p_base(new_p_base), number_types(new_p_base->numberTypes()), index_table(new_p_base->numberTypes()+1), elements(0)
	: p_base(new_p_base), number_types(new_p_base->numberTypes())
{
	int number_elements = 0;
	for (int j=0; j< number_types; ++j) {
		number_elements += p_base->numberStringsOfType(j);
	}
// 	elements.resize(number_elements);
	elements = new number[number_elements];
	index_table = new int[new_p_base->numberTypes()+1];
	index_table[0]=0;
	// look out! the last pointer points one past the allocated space. // useful to get size info.
	for (int i = 1; i <= number_types; ++i)
		index_table[i] = index_table[i-1] + p_base->numberStringsOfType(i-1);
	// zero out.
	for (int i = 0; i < number_elements; ++i)
		elements[i] = initial_value;
}


// copy constructor
template<class number>
Strip<number>::Strip (const Strip& original): p_base(original.p_base), number_types(original.numberTypes())
{
	index_table = new int[number_types+1];
	elements = new number[original.numberElements()];
	for (int i = 0; i <= number_types; ++i)
		index_table[i] = original.index_table[i];
	for (int i = 0; i < original.numberElements(); ++i)
		elements[i] = original.elements[i];
}



// assignment
template<class number>
Strip<number>& Strip<number>::operator= (const Strip<number>& original)
{
	if (numberElements() != original.numberElements()) throw Exception("Strip::operator=", exc_BadAssignASDEFINEDINTHISFILE, "(number elements)");
	if (number_types != original.numberTypes()) throw Exception("Strip::operator=", exc_BadAssignASDEFINEDINTHISFILE, "(number types)");

	if (this != &original) {
		for (int i = 0; i <= number_types; ++i) {
			if (index_table[i] != original.index_table[i])
				throw Exception("Strip::operator=", exc_BadAssignASDEFINEDINTHISFILE, "(index table)");
		}
		// copy elements
		for (int i = 0; i < original.numberElements(); ++i)
			elements[i] = original.elements[i];
	}

}


// destructor
template<class number>
Strip<number>::~Strip() {
	delete[] elements;
	delete[] index_table;
}



// outputting a Strip
template<class number>
ostream& operator<< (ostream& stream, const Strip<number>& out_to_be_put) {

	if (out_to_be_put.numberTypes() > 0)
		for (int i=0; i< out_to_be_put.numberTypes(); ++i) {
			for (int j=0; j < out_to_be_put.numberStringsOfType(i); ++j) {
 				stream << out_to_be_put(i, j) << SEP;
			}
			stream << ";" << SEP;
		}
	return stream;
}

// sorting a strip
template<class number, class Compare>
void sort(Strip<number>& strip_to_sort, Compare compare)
{
	for (int i=0; i<strip_to_sort.numberTypes(); ++i)
		sort (strip_to_sort.elements + strip_to_sort.index_table[i], strip_to_sort.elements+strip_to_sort.index_table[i+1], compare);
}




/*
/* FullStrip
 */

/*
template<class number>
FullStrip<number>::FullStrip (Base* new_p_base, const number initial_value) : p_base(new_p_base), dimension(new_p_base->numberTypes())
{
	// reserve space for the pointer table
	contents = new number*[dimension+1];

	int number_elements = 0;
	for (int j=0; j< dimension; ++j) {
		number_elements += p_base->numberStringsOfType(j) * p_base->chain.stringLength(j);
	}
	// reserve space for the elements
	contents[0] = new number[number_elements];

	// look out! the last pointer points one past the allocated space. // useful to get size info.
	for (int i = 1; i <= dimension; ++i) contents[i] = contents[i-1] + p_base->numberStringsOfType(i-1);

	// zero out.
	for (int i = 0; i < number_elements; ++i) *(contents[0] + i) = initial_value;
}


// copy constructor
template<class number>
FullStrip<number>::FullStrip (const FullStrip& original): p_base(original.p_base), dimension(original.numberTypes())
{
	contents= new number* [dimension+1];
	contents[0] = new number[original.numberElements()];
	// copy elements.
	for (int i = 0; i < original.numberElements(); ++i) *(contents[0] + i) = *(original[0] + i);
	// 'copy' pointer table, but with pointers to the new array
	for (int i = 1; i <= dimension; ++i) contents[i] = contents[0] + (original[i] - original[0]);
}

// assignment
template<class number>
FullStrip<number>& FullStrip<number>::operator= (const FullStrip<number>& original)
{
	if (this != &original) {
		// copy elements
		for (int i = 0; i < original.numberElements(); ++i) {
			*(contents[0] + i) = *(original.contents[0] + i);
		}
	}
}

// destructor
template<class number>
FullStrip<number>::~FullStrip()
{
	if (numberElements()) delete[] contents[0];
	delete[] contents;
}

*/


template<class number>
FullStrip<number>::FullStrip (const Base* const new_p_base, const number initial_value)
	: p_base(new_p_base), number_types(new_p_base->numberTypes())
{
	int number_elements = p_base->numberRoots();

	elements = new number[number_elements];
	index_table = new int[number_types+1];

	index_table[0]=0;
	// look out! the last pointer points one past the allocated space. // useful to get size info.
	for (int i = 1; i <= number_types; ++i)
		index_table[i] = index_table[i-1] + p_base->numberStringsOfType(i-1) * p_base->p_chain->stringLength(i-1);

	// zero out.
	for (int i = 0; i < number_elements; ++i)
		elements[i] = initial_value;
}


// copy constructor
template<class number>
FullStrip<number>::FullStrip (const FullStrip& original): p_base(original.p_base), number_types(original.numberTypes())
{
	elements = new number[original.numberElements()];
	index_table = new int[number_types+1];

	// copy elements.
	for (int i = 0; i < original.numberElements(); ++i)
		elements[i] = original.elements[i];
	// copy index table
	for (int i = 0; i <= number_types; ++i)
		index_table[i] =  original.index_table[i];
}

// assignment
template<class number>
FullStrip<number>& FullStrip<number>::operator= (const FullStrip<number>& original)
{
	if (this != &original) {
		// check index table
		for (int i = 0; i <= number_types; ++i) {
			if (index_table[i] != original.index_table[i])
				throw Exception("FullStrip::operator=", exc_BadAssignASDEFINEDINTHISFILE, " (index table)");
		}
		// copy elements
		for (int i = 0; i < original.numberElements(); ++i) {
			elements[i] = original.elements[i];
		}
	}
}


// destructor
template<class number>
FullStrip<number>::~FullStrip()
{
	delete[] elements;
	delete[] index_table;
}




#else

/** ********************* ********************** ************************* *********************** **/
// (ifdef STALINIST_RANGE_CHECK)
/* Failsafe versions:
 */



template<class number>
Strip<number>::Strip (Base* new_p_base, const number initial_value)
	: p_base(new_p_base), number_types(new_p_base->numberTypes()), index_table(new_p_base->numberTypes()+1), elements(0)
{
	try {
		int number_elements = 0;
		for (int j=0; j< number_types; ++j) {
			number_elements += p_base->numberStringsOfType(j);
		}
		elements.resize(number_elements);

		index_table[0]=0;
		// look out! the last pointer points one past the allocated space. // useful to get size info.
		for (int i = 1; i <= number_types; ++i)
			index_table.at(i) = index_table.at(i-1) + p_base->numberStringsOfType(i-1);
		// zero out.
		for (int i = 0; i < number_elements; ++i)
			elements.at(i) = initial_value;
	}
	catch (std::out_of_range) {
		throw Exception("Strip constructor", exc_OutOfRangeASDEFINEDINTHISFILE);
	}
}


// copy constructor
template<class number>
Strip<number>::Strip (const Strip& original): p_base(original.p_base), number_types(original.numberTypes()), index_table(original.index_table), elements(original.elements)
{
}



// assignment
template<class number>
Strip<number>& Strip<number>::operator= (const Strip<number>& original)
{
	if (numberElements() != original.numberElements()) throw Exception("Strip::operator=", exc_BadAssignASDEFINEDINTHISFILE, "(number elements)");
	if (number_types != original.numberTypes()) throw Exception("Strip::operator=", exc_BadAssignASDEFINEDINTHISFILE, "(number types)");
	try {
		if (this != &original) {
			// copy elements
			elements = original.elements;
// 			for (int i = 0; i < original.numberElements(); ++i) {
// 				elements.at(i) = original.elements.at(i);
// 			}
			if (index_table.size() != number_types+1) throw Exception("Strip::operator=", exc_OutOfRangeASDEFINEDINTHISFILE, "(index table)");
			for (int i = 0; i <= number_types; ++i) {
				if (index_table.at(i) != original.index_table.at(i))
					throw Exception("Strip::operator=", exc_BadAssignASDEFINEDINTHISFILE, "(index table)");
			}
		}
	}
	catch (std::out_of_range) {
		throw Exception("Strip::operator=", exc_OutOfRangeASDEFINEDINTHISFILE, "(index table)");
	}
}


// destructor
template<class number>
Strip<number>::~Strip() { }



// outputting a Strip
template<class number>
ostream& operator<< (ostream& stream, const Strip<number>& out_to_be_put) {

	if (out_to_be_put.numberTypes() > 0)
		for (int i=0; i< out_to_be_put.numberTypes(); ++i) {
			for (int j=0; j < out_to_be_put.numberStringsOfType(i); ++j) {
 				stream << out_to_be_put(i, j) << SEP;
			}
			stream << ";" << SEP;
		}
	return stream;
}

// sorting a strip
template<class number, class Compare>
void sort(Strip<number>& strip_to_sort, Compare compare)
{
	for (int i=0; i<strip_to_sort.numberTypes(); ++i)
		sort (strip_to_sort.elements.begin()+ strip_to_sort.index_table.at(i), strip_to_sort.elements.begin()+strip_to_sort.index_table.at(i+1), compare);
}




/// fullstrip

template<class number>
FullStrip<number>::FullStrip (Base* new_p_base, const number initial_value)
	: p_base(new_p_base), number_types(new_p_base->numberTypes()), index_table(new_p_base->numberTypes()+1), elements(0)
{
	try {
		int number_elements = 0;
		for (int j=0; j< number_types; ++j) {
			number_elements += p_base->numberStringsOfType(j) * p_base->chain.stringLength(j);
		}
		elements.resize(number_elements);

		index_table[0]=0;
		// look out! the last pointer points one past the allocated space. // useful to get size info.
		for (int i = 1; i <= number_types; ++i)
			index_table.at(i) = index_table.at(i-1) + p_base->numberStringsOfType(i-1)* p_base->chain.stringLength(i-1);

		// zero out.
		for (int i = 0; i < number_elements; ++i)
			elements.at(i) = initial_value;
	}
	catch (std::out_of_range) {
		throw Exception("FullStrip constructor", exc_OutOfRangeASDEFINEDINTHISFILE);
	}
}


// copy constructor
template<class number>
FullStrip<number>::FullStrip (const FullStrip& original): p_base(original.p_base), number_types(original.numberTypes())
{
	try {
		// copy elements.
		for (int i = 0; i < original.numberElements(); ++i)
			elements.at(i) = original.elements.at(i);
		// 'copy' pointer table, but with pointers to the new array
		for (int i = 0; i <= number_types; ++i)
			index_table.at(i) =  original.index_table.at(i);
	}
	catch (std::out_of_range) {
		throw Exception("FullStrip copy", exc_OutOfRangeASDEFINEDINTHISFILE);
	}
}

// assignment
template<class number>
FullStrip<number>& FullStrip<number>::operator= (const FullStrip<number>& original)
{
	if (numberElements() != original.numberElements()) throw Exception("FullStrip::operator=", exc_OutOfRangeASDEFINEDINTHISFILE, "(number elements)");
	if (number_types != original.numberTypes()) throw Exception("FullStrip::operator=", exc_OutOfRangeASDEFINEDINTHISFILE, "(number types)");
	try {
		if (this != &original) {
			// copy elements
			for (int i = 0; i < original.numberElements(); ++i) {
				elements.at(i) = original.elements.at(i);
			}
			for (int i = 0; i <= number_types; ++i) {
				if (index_table.at(i) != original.index_table.at(i))
					throw Exception("FullStrip::operator=", exc_OutOfRangeASDEFINEDINTHISFILE, " (index table)");
			}
		}
	}
	catch (std::out_of_range) {
		throw Exception("FullStrip::operator=", exc_OutOfRangeASDEFINEDINTHISFILE);
	}
}


// destructor
template<class number>
FullStrip<number>::~FullStrip() { }

#endif


#endif
