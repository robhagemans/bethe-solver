#include "square.h"

const string exc_IncompatibleSizes = "incompatible matrix sizes";
const string exc_SizeOne = "cannot chop a row and column from a size one matrix";


// Returns the ln of absolute value or norm of determinant of matrix a, through LU decomposition
//

// note: these are DESTRUCTIVE: you don't get your matrix back! (because of memory optimisation)
double lndetDestroy (Square<double>& mat)
{
	vector<int> index(mat.size());
	int sign_permutation;

	decomposeLU (mat, index, sign_permutation);
	// log( 1.0*sign_permutation ); only possible with complex numbers
	// but we're really doing log norm det.
	double log_det = 0.0;
	for (int j = 0; j < mat.size(); j++)
		// I've added an abs here to avoid taking the log of a negative number
		// if we need the sign, get it separately.
		log_det += log(fabs(mat[j][j]));

	return log_det;
}
double lndetDestroy (Square< complex<double> >& mat)
{
	vector<int> index(mat.size());
	int sign_permutation;

	decomposeLU (mat, index, sign_permutation);
	// log( 1.0*sign_permutation ); only possible with complex numbers
	// but we're really doing log norm det.
	double log_det = 0.0;
	for (int j = 0; j < mat.size(); j++)
		// I've added an abs here to avoid taking the log of a negative number
		// if we need the sign, get it separately.
		log_det += log(abs(mat[j][j]));

	return log_det;
}
