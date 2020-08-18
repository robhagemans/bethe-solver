Dynamics of Heisenberg Spin chains using the Bethe Ansatz
=========================================================

This repository contains code to find high-precision numerical solutions to the Bethe Ansatz equations for the Heisenberg spin chains in its three regimes:
- isotropic
- gapless aisotropic
- gapped

Further code is included to calculate form factors which can be used to find approximate solutions to the dynamical correlation functions for these models.

For the isotropic regime, additional code is provided to find exact solutions in cases where the Takahashi string hypothesis breaks down; in particular, deviated strings, extra real solutions and singular symmetric states.

This is research code; it is provided as-is, it is no longer maintained and parts of it may not work. It is shared here under a free licence in the hope that it may be useful for future students or researchers of these systems.

The algorithms underlying the code in this repository were developed and exposed in my doctoral thesis:
- R L Hagemans, [_Dynamics of Heisenberg spin chains_](https://pure.uva.nl/ws/files/4416871/52466_hagemans_thesis.pdf) (Amsterdam, 2007)

The techniques are also described in detail in:
- R Hagemans, J-S Caux, _Deformed strings in the Heisenberg model_, J Phys A: Math Theor **40** 14605 (2007)
- R Hagemans, J-S Caux, J M Maillet, _How to calculate correlation functions of Heisenberg chains_, in: _Lectures on the physics of highly correlated electron systems X_, eds. A Avella, F Mancini, AIP Conf Proc **846** 245 (2006)

This code and/or these techniques were also used, to greater or smaller extent, in:
- J-S Caux, R Hagemans, _The four-spinon dynamical structure factor of the Heisenberg chain_, J Stat Mech: Theor Exp (2006) P12013
- R G Pereira, J Sirker, J-S Caux, R Hagemans, J M Maillet, S R White, I Affleck, _The dynamical spin structure factor for the anisotropic spin-1/2 Heisenberg chain_, Phys Rev Lett **96** 257202 (2006)
- J-S Caux, R Hagemans, J M Maillet, _Computation of dynamical correlation functions of Heisenberg chains: the gapless anisotropic regime_, J Stat Mech: Theor Exp (2005) P09003
- R G Pereira, J Sirker,J-S Caux, R Hagemans, J M Maillet, S R White, I Affleck,
_Dynamical structure factor at small q for the XXZ spin-1/2 chain_, J Stat Mech: Theor Exp (2007) P08022
- J Sirker, R G Pereira, J-S Caux, R. Hagemans, J M Maillet, S R White, I Affleck, _Boson decay and the dynamical structure factor for the XXZ chain at finite magnetic field_, Phys B: Cond Mat **403** 1520 (2008)
