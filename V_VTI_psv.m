function [V]=V_VTI_psv(C,n) %11,13,33,55,66
pm=[1;-1];
n1=n(1);

C13p55=C(1,3)+C(5,5);
C33p55=C(3,3)+C(5,5);
C11m33=C(1,1)-C(3,3);
C33m55=C(3,3)-C(5,5);
C135  =C(1,1)+C(3,3)-2*C(5,5);

%compute phase velocities
D1 =2*(2*C13p55^2-C33m55*C135);  
D2 =C135^2-4*C13p55^2;
D  =sqrt(C33m55^2+D1*n1^2+D2*n1^4);

V=zeros(2,1); 
for w=1:2
    V(w)=1/sqrt(2)*sqrt(C33p55+C11m33*n1^2+pm(w)*D);   %phase velocity
end
return
