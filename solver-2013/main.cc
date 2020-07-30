#include "exception.h"

string Command;
vector< stringstream* > Arguments;
stringstream Command_Line;

int main(int argc, char* argv[])
{
	Command = argv[0];
	for (int argn=1; argn < argc; ++argn) {
		Command_Line << argv[argn] << "\t";
		Arguments.push_back(new stringstream(argv[argn]));
	}
	
	try {	
		return run();
	}
	catch (Exception error) {
		cerr << "Error: " << error<< endl;
		return 1;
	}
	catch (const char* error) {
		cerr << "Error: " << error<< endl;
		return 1;
	}
	catch (char* error) {
		cerr << "Error: " << error<< endl;
		return 1;
	}
	catch (int error) {
		cerr << "Error: " << error<< endl;
		return 1;
	}
	catch (const int error) {
		cerr << "Error: " << error<< endl;
		return 1;
	}
	catch (std::bad_alloc) {
		cerr << "Error: badalloc"<<endl;
		return 1;
	}
	catch (std::bad_exception) {
		cerr << "Error: badexception"<<endl;
		return 1;
	}
	catch (std::out_of_range) {
		cerr << "Error: index out of range"<<endl;
		return 1;
	}
	catch (...) {
		cerr << "Unknown error. Sorry 'bout that. "<<endl;
		return 1;
	}
	
	for (int i=0; i < Arguments.size(); ++i) delete Arguments[i];
}
