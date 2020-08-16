function [J0, J1, J2, J3, J4, J5, J6, J7] = quantum(N, l0,l1,l2,l3,l4,l5, l6, l7)
    J0 = (1/pi) * (N*atan(2*l0) - atan(l0-l1) - atan(l0-l2) - atan(l0-l3) - atan(l0-l4) - atan(l0-l5) - atan(l0-l6) - atan(l0-l7));    
    J1 = (1/pi) * (N*atan(2*l1) - atan(l1-l0) - atan(l1-l2) - atan(l1-l3) - atan(l1-l4) - atan(l1-l5) - atan(l1-l6) - atan(l1-l7));    
    J2 = (1/pi) * (N*atan(2*l2) - atan(l2-l1) - atan(l2-l0) - atan(l2-l3) - atan(l2-l4) - atan(l2-l5) - atan(l2-l6) - atan(l2-l7));    
    J3 = (1/pi) * (N*atan(2*l3) - atan(l3-l1) - atan(l3-l2) - atan(l3-l0) - atan(l3-l4) - atan(l3-l5) - atan(l3-l6) - atan(l3-l7));    
    J4 = (1/pi) * (N*atan(2*l4) - atan(l4-l1) - atan(l4-l2) - atan(l4-l3) - atan(l4-l0) - atan(l4-l5) - atan(l4-l6) - atan(l4-l7));    
    J5 = (1/pi) * (N*atan(2*l5) - atan(l5-l1) - atan(l5-l2) - atan(l5-l3) - atan(l5-l4) - atan(l5-l0) - atan(l5-l6) - atan(l5-l7));    
    J6 = (1/pi) * (N*atan(2*l6) - atan(l6-l1) - atan(l6-l2) - atan(l6-l3) - atan(l6-l4) - atan(l6-l5) - atan(l6-l0) - atan(l6-l7));    
    J7 = (1/pi) * (N*atan(2*l7) - atan(l7-l1) - atan(l7-l2) - atan(l7-l3) - atan(l7-l4) - atan(l7-l5) - atan(l7-l6) - atan(l7-l0));    
endfunction

#function [J0, J1, J2, J3] = quantum(N, l0,l1,l2,l3)
#    J0 = (1/pi) * (N*atan(2*l0) - atan(l0-l1) -atan(l0-l2) - atan(l0-l3)  );    
#    J1 = (1/pi) * (N*atan(2*l1) - atan(l1-l0) -atan(l1-l2) - atan(l1-l3)  );    
#    J2 = (1/pi) * (N*atan(2*l2) - atan(l2-l1) -atan(l2-l0) - atan(l2-l3)  );    
#    J3 = (1/pi) * (N*atan(2*l3) - atan(l3-l1) -atan(l3-l2) - atan(l3-l0)  );    
#        
#endfunction
