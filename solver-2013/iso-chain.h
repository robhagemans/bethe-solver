#ifndef ISO_CHAIN_H
#define ISO_CHAIN_H

#include "chain.h"

// maximum number of string types to consider
/// NOTE: number lower than 4 screw things up.
#define CUTOFF_TYPES 4

	
	/** XXX Chain: **/
	
	class XXX_Chain : public ::Chain {
	public:
		// create 
		XXX_Chain(const int chain_length, const int cutoff_types=CUTOFF_TYPES);
		// copy
		XXX_Chain(const Chain& original): Chain(original) { if (original.delta() != 1.0) throw Exception("XXXChain copy", exc_CopyMismatch); }; 
		// destructor
		virtual ~XXX_Chain(void) {};
		// clone
		virtual inline Chain* clone(void) const { return new XXX_Chain(*this); };
		// returns anisotropy delta
		virtual inline REAL delta(void) const { return 1.0;};
	};


#endif
