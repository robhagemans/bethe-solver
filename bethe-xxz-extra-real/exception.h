#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
using namespace std;



extern string Command;
extern vector< stringstream* > Arguments;
extern stringstream Command_Line;

int run(void);


class Exception {
public:
	string location;
	string error;
	string remarks;

	Exception(string here, string message);
	Exception(string here, string message, string sowhat);
};

ostream& operator<< (ostream& stream, const Exception& out_to_be_put);

#endif
