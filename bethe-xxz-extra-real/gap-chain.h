#ifndef GAP_CHAIN_H
#define GAP_CHAIN_H

#include "exception.h"
#include "bethe.h"
#include "chain.h"


	/** Gapped Chain: **/

	class Gap_Chain : public ::Chain {
	public:
		// construct
		Gap_Chain(const REAL delta_anisotropy, const int chain_length, const int cutoff_types=100);
		// copy
		Gap_Chain(const Chain& original): Chain(original) { if (original.delta() <= 1.0) throw Exception("Gap_chain copy", exc_CopyMismatch); };
		// destructor
		virtual ~Gap_Chain(void);
		// clone
		virtual inline Chain* clone(void) const { return new Gap_Chain(*this); };
		// returns anisotropy delta
		virtual inline REAL delta(void) const;
	};

#endif
