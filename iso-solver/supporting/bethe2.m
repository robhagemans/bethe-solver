function [l0p,epp,l2p,l3p] = bethe2(l0,ep,l2,l3)
    N=12; J0=4.5; J1=3.5; J2=4.5; J3=5.5;
    l1=real(l0+ep);
    #epp=real(tan( N*atan(2*l1)-pi*J1-atan(l1-l2)-atan(l1-l3) ));
    epp=real(tan( 0.5*( N*atan(2*l1)-N*atan(2*l0)  -pi*(J1-J0) -atan(l1-l2)-atan(l1-l3)+atan(l0-l2)+atan(l0-l3) )));
    
    l0p= real(0.5*tan( (1.0/N) * (pi*J0 + atan(-ep) + atan(l0-l2) + atan(l0-l3)) ));
    l2p= 0.5*tan( (1.0/N) * (pi*J2 + atan(l2-l0) + atan(l2-l1) + atan(l2-l3)) );
    l3p= conj(l2p)
 endfunction
