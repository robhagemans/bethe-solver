#include "configuration.h"

IsoConfiguration::IsoConfiguration(const int number_roots)
	: number_strings_(1, number_roots)
{	}

IsoConfiguration::IsoConfiguration(const vector<int>& base_vec)
	: number_strings_(base_vec)
{	}

IsoConfiguration::IsoConfiguration(const int* base_array, const int array_length)
	: number_strings_(base_array, base_array+array_length)
{	}


// add strings to an existing configuration. add negative numbers to remove strings.
// if number of strings ends up less than zero, sets to zero
IsoConfiguration& IsoConfiguration::add(const int string_type, const int number_strings_to_add)
{
    if (number_strings_.size() < string_type+1) {
        number_strings_.resize(string_type+1, 0);
    }
    number_strings_[string_type] += number_strings_to_add;
    if (number_strings_[string_type] < 0) {
        number_strings_[string_type] = 0;
    }

    return *this;
}

// set the number of strings of a given type to a number
// if number of strings ends up less than zero, sets to zero
IsoConfiguration& IsoConfiguration::set(const int string_type, const int number_strings)
{
    if (number_strings_.size() < string_type+1) {
        number_strings_.resize(string_type+1, 0);
    }
    number_strings_[string_type] = number_strings;
    if (number_strings_[string_type] < 0) {
        number_strings_[string_type] = 0;
    }

    return *this;
}



inline int IsoConfiguration::stringLength(const int type) const
{
	if (type <0) return 0;
	if (type >= number_strings_.size()) return 0;
	return type+1;
}

inline vector<int> IsoConfiguration::stringLengths() const
{
	vector<int> lengths (number_strings_.size());
	for (int i=0; i< lengths.size(); ++i)
		lengths[i] = i+1;
	return lengths;
}

inline int IsoConfiguration::numberStringsOfType(const int type)  const
{
	if (type < 0) return 0;
	if (type >= number_strings_.size()) return 0;
	return number_strings_[type];
}

inline int IsoConfiguration::numberTypes()  const
{
	return number_strings_.size();
}

inline int IsoConfiguration::numberRoots() const
{
	int number_roots = 0;
	for (int i=0; i<number_strings_.size(); ++i)
		number_roots += number_strings_[i]*stringLength(i);
	return number_roots;
}


/*
inline int IsoConfiguration::numberStringsOfLength(const int length) const
{
	return numberStringsOfType(length-1);
}
*/


/*

// add string of given length (length, not type!)
inline bool IsoConfiguration::insertStringOfLength(const int length)
{
	if (length < 1) return false;
	if (length > number_strings_.size()) number_strings_.resize(length);
	++number_strings_[length-1];
	return true;
}

// remove string of given length (length, not type!); false if no string present
inline bool IsoConfiguration::deleteStringOfLength(const int length)
{
	if (length < 1) return false;
	if (length > number_strings_.size()) return false;
	if (number_strings_[length-1] <1) return false;
	--number_strings_[length-1];
	return true;
}
*/
