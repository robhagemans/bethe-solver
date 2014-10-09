#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <vector>

using namespace std;


/** configuration ABC **/


class Configuration {
public:
	// length of string of type number, zero if not defined
	virtual int stringLength(const int type) const = 0;
	// length of string of type number
	virtual vector<int> stringLengths() const = 0;
	// number strings of a given type number
	virtual int numberStringsOfType(const int type) const = 0;
	// number of string types in configuration
	virtual int numberTypes() const = 0;
	// total number of roots 
	virtual int numberRoots() const = 0;

};


/** configuration **/

class IsoConfiguration: public Configuration {

protected:
	vector<int> number_strings_;
public:
	// default costructor: trivial configuration
	IsoConfiguration(const int number_magnons=0);
	// construct from given vector/array of number of strings for each string type 
	IsoConfiguration(const vector<int>& base_vec);
	IsoConfiguration(const int* base_array, const int array_length);
	~IsoConfiguration() {};	

    // add strings to an existing configuration. add negative numbers to remove strings.
    IsoConfiguration& add(const int string_type, const int number_strings_to_add);
    // set the number of strings of a given type to a number
    IsoConfiguration& set(const int string_type, const int number_strings);

	// number strings of type number
	int numberStringsOfType(const int type) const;
	// number of string types in configuration
	int numberTypes() const;
	// total number of roots 
	int numberRoots() const;

	// length of string of type number
	int stringLength(const int type) const;
	vector<int> stringLengths() const;
	
};




#endif
