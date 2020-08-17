#include "exception.h"


Exception::Exception(string here, string message)
	: where_(here), what_(message), why_("")
{	}

Exception::Exception (string here, string message, string so_what)
	: where_(here), what_(message), why_(so_what)
{	}


string Exception::why() const
{
	return why_;
}

string Exception::what() const
{
	return what_;
}

string Exception::where() const
{
	return where_;
}



ostream& operator<< (ostream& stream, const std::exception& out_to_be_put)
{
	stream<< out_to_be_put.what();
	return stream;
}

ostream& operator<< (ostream& stream, const Exception& out_to_be_put)
{
	stream<< out_to_be_put.where() <<": "<<out_to_be_put.what()<<"; "<<out_to_be_put.why();
	return stream;
}
