function [l0p,l1p,l2p,l3p] = bethe(N, J0,J1,J2,J3, l0,l1,l2,l3)
    l0p= 0.5*tan( (1.0/N) * (pi*J0 + atan(l0-l1) + atan(l0-l2) + atan(l0-l3)) );
    l1p= 0.5*tan( (1.0/N) * (pi*J1 + atan(l1-l0) + atan(l1-l2) + atan(l1-l3)) );
    l2p= 0.5*tan( (1.0/N) * (pi*J2 + atan(l2-l0) + atan(l2-l1) + atan(l2-l3)) );
    l3p= 0.5*tan( (1.0/N) * (pi*J3 + atan(l3-l0) + atan(l3-l1) + atan(l3-l2)) );
endfunction

