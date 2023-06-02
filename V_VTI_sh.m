function [V]=V_VTI_sh(C,n) 
n1=n(1);n2=n(2);
V=sqrt(C(6,6)*n1^2+C(5,5)*n2^2);
return 
