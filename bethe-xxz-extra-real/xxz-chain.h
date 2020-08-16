#ifndef XXZ_CHAIN_H
#define XXZ_CHAIN_H

#include "exception.h"
#include "bethe.h"
#include "chain.h"

#define CONT_FRAC_DEFAULT_PRECISION 1e-6
#define CONT_FRAC_CUTOFF 100


extern const char* exc_Rational;

	/** XXZ Chain: **/

	class XXZ_Chain : public ::Chain {
	public:
		// constructor
		XXZ_Chain(const REAL delta_anisotropy, const int chain_length, const REAL precision=CONT_FRAC_DEFAULT_PRECISION, const int cutoff_fraction=CONT_FRAC_CUTOFF);
		// copy constructor
		XXZ_Chain(const Chain& original): Chain(original) { if (original.delta() >= 1.0 || original.delta() < 0.0) throw Exception("XXZChain copy", exc_CopyMismatch); };
		// destructor
		virtual ~XXZ_Chain(void) {};
		// clone
		virtual inline Chain* clone(void) const { return new XXZ_Chain(*this); };
		// returns anisotropy delta
		virtual inline REAL delta(void) const {	return cos(anisotropy);	};
	};

#endif
