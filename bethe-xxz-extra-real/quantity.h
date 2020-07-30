#ifndef QUANTITY_H
#define QUANTITY_H


/** magnetic fields **/
// a table of magnetisations and magnetic fields (for finite-temperature calculations)
vector<REAL> magneticFields (const Chain* const p_chain);
// the magnetic field for a given chain and magnetisation // this is needed to calculate the energy for S-minus matrix elements
REAL magneticField (const Chain* const p_chain, const int number_down);



/** class Quantity: form factor, max sum and other form-factor dependent functions for various chains **/
class Quantity {
protected:
	State* left_state; // ground state
	State* right_state;
	mutable REAL form_factor;
	mutable REAL rotcaf_mrof;
public:
	// default constructor: empty quantity
	inline Quantity (void) : left_state(0), right_state(0), form_factor(NOT_CALCULATED), rotcaf_mrof(NOT_CALCULATED){};
	// set only left state
	inline Quantity (State* p_left) : left_state(p_left), right_state(0), form_factor(NOT_CALCULATED), rotcaf_mrof(NOT_CALCULATED){};
	// set both left and right states
	inline Quantity (State* p_left, State* p_right) : left_state(p_left), right_state(p_right), form_factor(NOT_CALCULATED), rotcaf_mrof(NOT_CALCULATED){};
	
	// (re-)set left state
	inline void setLeftState (State* p_left) { left_state = p_left; form_factor = NOT_CALCULATED; rotcaf_mrof = NOT_CALCULATED;};
	// (re-)set right state
	inline void setRightState (State* p_right) { right_state = p_right; form_factor = NOT_CALCULATED; rotcaf_mrof = NOT_CALCULATED;};
	
	// access right state
	inline State* pRightState (void) const { return right_state; };
	// access left state
	inline State* pLeftState (void) const { return left_state; };
	
	// form factor sq ( < left | op | right > )
	inline REAL formFactor (void) const	{ return form_factor==NOT_CALCULATED ? calculateFormFactor() : form_factor; };
	// form factor sq ( < right | op | left > )
	inline REAL rotcafMrof (void) const { return rotcaf_mrof==NOT_CALCULATED ? rotcafMrofEtaluclac() : rotcaf_mrof; };
	
	// sum rule with given chain length and number of down spins
	virtual REAL maxSum (const int number_sites, const int left_number_down) const = 0; 
	// sum rule with current left state
	inline REAL maxSum (void) const { return maxSum (left_state->p_base->p_chain->length(), left_state->p_base->numberDown()); };
	
	// sum rule excluding infinite rapidities (xxx)
	virtual REAL maxSumHighestWeight (const int number_sites, const int left_number_down) const = 0; 
	// sum rule excluding infinite rapidities with current left state
	inline REAL maxSumHighestWeight (void) const { return maxSumHighestWeight (left_state->p_base->p_chain->length(), left_state->p_base->numberDown()); };
	
	// number of magnons in left state
	inline int leftNumberDown (void) const { return left_state->p_base->numberDown(); };
	// number of magnons in right state that gives nonzero contribution if the left state has a given number of magnons
	virtual int rightNumberDown (int left_number_down) const = 0; 
	// number of magnons in right state that gives nonzero contribution with current left state
	inline int rightNumberDown (void) const { return rightNumberDown(left_state->p_base->numberDown()); };
	
	// energy shift due to magnetic field
	virtual REAL energyShift (void) const = 0;
	// energy with respect to left (ground) state, including field contribution
	inline REAL energy (void) const { return right_state->energy() - left_state->energy() + energyShift(); } 
	
	// momentum index with respect to left state	
	inline int mode (void) const { return modulo(right_state->mode() - left_state->mode(), left_state->p_chain->length()); }
	// momentum with respect to left state	
	inline REAL momentum (void) const { return 2.0*PI*mode()/(1.0*left_state->p_chain->length()); }
	
	// name of quantity, for file names
	virtual string name (void) const = 0;
	// name of related structure factor, for file names
	virtual string sfName (void) const = 0;
	
	// maximum number of infinite rapidities for which form factor is nonzero
	virtual int maxInfinite (void) const = 0;
	// whether it makes sense to deviate (decause we can't use deviated state as right state)
	virtual bool deviable (void) const = 0;
	// base name for files relating to a calculation of this quantity
	inline const char* fileName (const string extension="") const 
	{ 	return (name() + ::name(left_state->p_chain, left_state->p_base->numberDown()) +"."+ extension).c_str(); };
	
protected:
	virtual REAL calculateFormFactor (void) const = 0;
	virtual REAL rotcafMrofEtaluclac (void) const = 0;
};


/** dummy quantity: calculate no form factors **/
class NullQuantity : public Quantity {
public:
	inline NullQuantity (void) : Quantity() {};
	inline NullQuantity (State* const p_left) : Quantity(p_left) {};
	inline NullQuantity (State* const p_left, State* const p_right) : Quantity(p_left,p_right) {};
	
	// sum rules are meaningless here
	inline virtual REAL maxSum (const int number_sites, const int left_number_down) const {	return 1.0; } ;
	inline virtual REAL maxSumHighestWeight (const int number_sites, const int left_number_down) const { return 1.0; };
	
	// both sides use same number of magnons
	inline virtual int rightNumberDown (const int left_number_down) const { return left_number_down; };
	inline virtual REAL energyShift (void) const { return 0.0; };
	
	// infinite rapidities: whatever.
	virtual inline int maxInfinite(void) const { return (left_state->p_chain->delta() == 1.0)?leftNumberDown():0; };
	// does it make sense to deviate? yes, if XXX
	virtual inline bool deviable(void) const { return (right_state->p_chain->delta()==1.0); };
	
	inline string name (void) const { return ""; } ;
	inline string sfName (void) const { return ""; };
protected:
	inline REAL calculateFormFactor (void) const { return form_factor = 0.0; } ; 
	inline REAL rotcafMrofEtaluclac (void) const { return rotcaf_mrof = 0.0; } ;	

};

/** longitudinal form factor **/
class Longitudinal : public Quantity {
public:
	static const string ff_name;
	static const string sf_name;
public:
	inline Longitudinal (void) : Quantity() {};
	inline Longitudinal (State* const p_left) : Quantity(p_left) {};
	inline Longitudinal (State* const p_left, State* const p_right) : Quantity(p_left,p_right) {};
	
	// full sum rule for Szz
	inline virtual REAL maxSum (const int number_sites, const int left_number_down) const 
	{	return 0.25*number_sites*(  1.0 -sq( 1.0-2.0*left_number_down/(1.0*number_sites) )  ); } ;
	// sum rule excluding infinite rapidities (xxx)
	inline virtual REAL maxSumHighestWeight (const int number_sites, const int left_number_down) const
	{	return 	0.25*number_sites*( 1.0 -sq(1.0-2.0*left_number_down/(1.0*number_sites)) ) - 1.0*left_number_down/REAL(number_sites-2.0*(left_number_down-1));	};
	// both sides need sme number of magnons
	inline virtual int rightNumberDown (const int left_number_down) const { return left_number_down; };
	// no energy shift due to field as both sides have the same number of magnons
	inline virtual REAL energyShift (void) const { return 0.0; };
	// infinite rapidities: yes, one, if we're on XXX.
	virtual inline int maxInfinite(void) const { return (left_state->p_chain->delta() == 1.0)?1:0; };
	// does it make sense to deviate?
	virtual inline bool deviable(void) const { return (right_state->p_chain->delta()==1.0) && (0==right_state->p_base->numberInfiniteRapidities()); };
	
	inline string name (void) const { return ff_name; } ;
	inline string sfName (void) const { return sf_name; };
	
protected:
	// NOTE that this may actually be calculated using the converse form factor!
	inline REAL calculateFormFactor (void) const { return form_factor = right_state->longitudinalFormFactor(*left_state); } ;
	inline REAL rotcafMrofEtaluclac (void) const { return rotcaf_mrof = left_state->longitudinalFormFactor(*right_state); } ;	
};


/** transverse (S-+) form factor **/
class Transverse : public Quantity {
protected:
	mutable REAL its_energy_shift;
public:
	static const string ff_name;
	static const string sf_name;
public:
	inline Transverse (void) : Quantity(), its_energy_shift(NOT_CALCULATED) {};
	inline Transverse (State* const  p_left) : Quantity(p_left), its_energy_shift(NOT_CALCULATED) {};
	inline Transverse (State* const p_left, State* const p_right) : Quantity(p_left, p_right), its_energy_shift(NOT_CALCULATED) {};
	inline Transverse (State* const p_left, State* const p_right, const REAL the_energy_shift) : Quantity(p_left, p_right), its_energy_shift(the_energy_shift) {};
	
	// full sum rule for S-+ 
	inline virtual REAL maxSum (const int number_sites, const int left_number_down) const { return left_number_down; };
	// infinite rapidities give zero for S-+
	inline virtual REAL maxSumHighestWeight (const int number_sites, const int left_number_down) const
	{	return Transverse::maxSum(number_sites,left_number_down); };
	// reference state has one less magnon
	inline virtual int rightNumberDown (const int left_number_down) const { return left_number_down-1; };
	// state and reference state have different magnetisations and are hence shifted in energy by the magnetic field
	inline virtual REAL energyShift (void) const
	{ 	return its_energy_shift==NOT_CALCULATED? its_energy_shift = -magneticField(left_state->p_base->p_chain, left_state->p_base->numberDown()) : its_energy_shift; };
	// no infinite rapidities in S- form factor
	virtual inline int maxInfinite(void) const { return 0; };
	// deviation makes no sense for S- as we always need the right state
	virtual inline bool deviable(void) const { return false; };

	inline string name (void) const { return ff_name; } ;
	inline string sfName (void) const { return sf_name; };
	
protected:
	inline REAL calculateFormFactor (void) const { return form_factor = right_state->transverseFormFactor(*left_state); } ;
	inline REAL rotcafMrofEtaluclac (void) const { return rotcaf_mrof = left_state->plusFormFactor(*right_state); } ;	
};



/** S+- form factor: implemented as S-+ the other way around **/
class SPlusMinus : public Quantity {
protected:
	mutable REAL its_energy_shift;
public:
	static const string ff_name;
	static const string sf_name;
public:
	inline SPlusMinus (void) : Quantity(), its_energy_shift(NOT_CALCULATED) {};
	inline SPlusMinus (State* const p_left) : Quantity(p_left), its_energy_shift(NOT_CALCULATED) {};
	inline SPlusMinus (State* const p_left, State* const p_right) : Quantity(p_left, p_right), its_energy_shift(NOT_CALCULATED) {};
	inline SPlusMinus (State* const p_left, State* const p_right, const REAL the_energy_shift) : Quantity(p_left, p_right), its_energy_shift(the_energy_shift) {};
	
	// full sum rule for S+- excluding the 'connected' part (i.e. the contribution from lower-weight states that have connected Szz form factor)
	// including connected part: {	return number_sites - left_number_down; };
	inline virtual REAL maxSum (const int number_sites, const int left_number_down) const 
 	{	return REAL(number_sites-left_number_down) - (1.0 - 2.0*left_number_down/REAL(number_sites));	};
	// sum rule excluding infinite rapidities (xxx)	(which automatically means excluding 'connected' stuff).
	inline virtual REAL maxSumHighestWeight (const int number_sites, const int left_number_down) const
	{	
		return 
	// full s+- sum rule
			REAL(number_sites - left_number_down) 
	// minus 1-inf contribution, including connected contributions. 
 			- 1.0/(1.0 - (2.0*left_number_down)/REAL(number_sites))
  			+ 4.0*left_number_down/REAL((number_sites-2*left_number_down)*(number_sites-2*(left_number_down-1)))
	// minus 2-inf contribution
			- 2.0*left_number_down/REAL((number_sites-2*left_number_down+2)*(number_sites-2*left_number_down+1));		
	};
	// reference state has one more magnon
	inline virtual int rightNumberDown (const int left_number_down) const { return left_number_down+1; };
	// state and reference state have different magnetisations and are hence shifted in energy by the magnetic field
	inline virtual REAL energyShift (void) const
	{ 	return its_energy_shift==NOT_CALCULATED? its_energy_shift = magneticField(left_state->p_base->p_chain, left_state->p_base->numberDown()) : its_energy_shift; };
	
	// infinite rapidities: yes, two, if we're on XXX.
	virtual inline int maxInfinite(void) const { return (left_state->p_chain->delta() == 1.0)?2:0; };
	// whether it makes sense to deviate
	virtual inline bool deviable(void) const { return (right_state->p_chain->delta()==1.0) && (right_state->p_base->numberInfiniteRapidities()<2); };

	inline string name (void) const { return ff_name; } ;
	inline string sfName (void) const { return sf_name; };
		
protected:
	inline REAL calculateFormFactor (void) const { return form_factor = right_state->plusFormFactor(*left_state); } ; 
	inline REAL rotcafMrofEtaluclac (void) const { return rotcaf_mrof = left_state->transverseFormFactor(*right_state); } ;	
};



/** generic constructors **/
// create from function index
inline Quantity* newQuantity(const int function_index, State* const p_left) 
{	
	if (function_index==1) return new Transverse(p_left); 
	else if (function_index==-1) return new SPlusMinus(p_left); 
	else if (function_index==0) return new Longitudinal(p_left);  
	else return new NullQuantity (p_left);
};
// create from name
inline Quantity* newQuantityFromName(const string quantity_name, State* const p_left) 
{ 	
	if (quantity_name==Longitudinal::ff_name) return new Longitudinal(p_left); 
	else if (quantity_name==Transverse::ff_name) return new Transverse(p_left); 
	else if (quantity_name==SPlusMinus::ff_name) return new SPlusMinus(p_left); 
	else return new NullQuantity (p_left); 
};


/** convert a (file) name into a Quantity **/
Quantity* newQuantityFromFileName (const string name, State* const p_ground_state);


/** generate list of bases **/
inline vector<BaseData> allNewBases (const Quantity* const p_quantity, const int max_string_length, const int uptoinc_number_particles, const int uptoinc_number_spinons)
{ return allNewBases (	p_quantity->pLeftState()->p_chain, p_quantity->rightNumberDown(), 
						max_string_length, uptoinc_number_particles, uptoinc_number_spinons,
						0	//p_quantity->maxInfinite() 	// zero if we do separate lower-weight scans.
					); 	};




#endif
