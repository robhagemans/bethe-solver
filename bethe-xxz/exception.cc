#include "exception.h"

const char* exc_IndexRange = "index out of range";



Exception::Exception(string thrown_here, string message)
{
	location = thrown_here;
	error = message;
	remarks = "";
}

Exception::Exception (string here, string message, string sowhat)
{
	location = here;
	error = message;
	remarks = sowhat;
}

ostream& operator<< (ostream& stream, const Exception& out_to_be_put) 
{
	stream<< out_to_be_put.location <<": "<<out_to_be_put.error<<" "<<out_to_be_put.remarks;
	return stream;
}
