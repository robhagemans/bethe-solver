#ifndef CHAIN_H
#define CHAIN_H

#include "bethe.h"

extern const char* exc_OddChain;
extern const char* exc_CopyMismatch;


class Chain {
	protected:
		// anisotropy *parameter* . This is zeta=cf(delta).
		REAL anisotropy; 
		int number_sites;
		int number_types;
	public:
// 	protected:
		int* string_length;
		int* string_parity;	
	public:	
		// empty constructor, let descendants set the fields (TODO: this is not quite optimal, of course...)
		Chain (void) {};
		// copy
		Chain (const Chain& original);
		// destructor
		virtual ~Chain(void);
		// make a new pointer and copy ourselves to it 
		virtual Chain* clone(void) const = 0;
		
		// same anisotropy and number sites? same chain.
		inline bool operator== (const Chain& another) const { return ( (number_sites==another.number_sites) && (anisotropy==another.anisotropy) );  };
		// not same chain?
		inline bool operator!= (const Chain& another) const { return !(*this==another); };
		
		// the physical anisotropy delta
		virtual REAL delta(void) const =0; 
		// anisotropy parameter zeta
		inline REAL zeta (void) const { return anisotropy; };
		// number of sites on the chain
	    inline int length (void) const { return number_sites; };
	    // number of recognised string types
		inline int numberTypes (void) const { return number_types; };
		// string length
		inline int stringLength (const int string_type) const { return string_length[string_type]; };
		// string parity
		inline int stringParity (const int string_type) { return string_parity[string_type]; };
		// is the parity positive?
		inline bool positiveParity (const int string_type) const { return (string_parity[string_type]==1); };
	};
	
#endif
