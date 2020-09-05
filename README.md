Dynamics of Heisenberg spin chains using the Bethe Ansatz
=========================================================

This repository contains code to find high-precision numerical solutions to the Bethe Ansatz equations for the Heisenberg spin chains in its three regimes:
- isotropic
- gapless aisotropic
- gapped

Further code is included to calculate form factors which can be used to find approximate solutions to the dynamical correlation functions for these models.

For the isotropic regime, additional code is provided to find exact solutions in cases where the Takahashi string hypothesis breaks down; in particular, deviated strings, extra real solutions and singular symmetric states.

This is research code. It is provided as-is and is no longer maintained. It may or may not work. In fact, as-is it does not work as you'll still need to add three functions from Numerical Recipes in C (`polint`, `ludcmp`, `lubksb`).

This code is shared here under a free licence in the hope that it may be useful for future students or researchers of these systems, but I am not holding my breath.

The algorithms underlying the code in this repository were developed and exposed in my doctoral thesis:
- R L Hagemans, [_Dynamics of Heisenberg spin chains_](https://pure.uva.nl/ws/files/4416871/52466_hagemans_thesis.pdf) (Amsterdam, 2007)

The techniques are also described in detail in:
- R Hagemans, J-S Caux,  
  _Deformed strings in the Heisenberg model_,  
  J Phys A: Math Theor **40** 14605 (2007), [arXiv:0707.2803](https://arxiv.org/abs/0707.2803)
- R Hagemans, J-S Caux, J M Maillet,  
  _How to calculate correlation functions of Heisenberg chains_,  
  in: A Avella, F Mancini (eds) _Lectures on the physics of highly correlated electron systems X_,  
  AIP Conf Proc **846** 245 (2006), [cond-mat/0611467](https://arxiv.org/abs/cond-mat/0611467)

This code and/or these techniques were also used, to greater or smaller extent, in:
- J-S Caux, R Hagemans, J M Maillet,  
  _Computation of dynamical correlation functions of Heisenberg chains: the gapless anisotropic regime_,  
  J Stat Mech: Theor Exp (2005) P09003, [cond-mat/0506698](https://arxiv.org/abs/cond-mat/0506698)
- J-S Caux, R Hagemans,  
  _The four-spinon dynamical structure factor of the Heisenberg chain_,  
  J Stat Mech: Theor Exp (2006) P12013, [cond-mat/0611319](https://arxiv.org/abs/cond-mat/0611319)
- R G Pereira, J Sirker, J-S Caux, R Hagemans, J M Maillet, S R White, I Affleck,  
  _The dynamical spin structure factor for the anisotropic spin-1/2 Heisenberg chain_,  
  Phys Rev Lett **96** 257202 (2006), [cond-mat/0603681](https://arxiv.org/abs/cond-mat/0603681)
- R G Pereira, J Sirker, J-S Caux, R Hagemans, J M Maillet, S R White, I Affleck,  
  _Dynamical structure factor at small q for the XXZ spin-1/2 chain_,  
  J Stat Mech: Theor Exp (2007) P08022, [arXiv:0706.4327](https://arxiv.org/abs/0706.4327)
- J Sirker, R G Pereira, J-S Caux, R Hagemans, J M Maillet, S R White, I Affleck,  
  _Boson decay and the dynamical structure factor for the XXZ chain at finite magnetic field_,  
  Phys B: Cond Mat **403** 1520 (2008), [arXiv:0705.3312](https://arxiv.org/abs/0705.3312)
