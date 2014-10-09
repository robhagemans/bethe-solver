#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
using namespace std;


class Exception {
protected:
	string what_;
	string where_;
	string why_;
public:
	virtual string what() const;
	virtual string why() const;
	virtual string where() const;

	Exception(const string here, const string message);
	Exception(const string here, const string message, string sowhat);
	
	~Exception() {};
};

ostream& operator<< (ostream& stream, const std::exception& out_to_be_put);
ostream& operator<< (ostream& stream, const Exception& out_to_be_put);

#endif
