%compute sin of Snell's phase angle from input sin phase angle and phase velocity of one layer
function [mu]=refx(w,u,V,C) %w=wavetype
if(u==0)          %vertical incidence, in case invhs is infinity
    mu=0;
else
    hs=u/V;       %horizontal slowness
    invhs=2/hs^2;
    a=(invhs)^2-2*C(1,1)*invhs+2*C(3,3)*invhs+4*C(1,3)^2-4*C(1,1)*C(3,3)+4*C(1,1)*C(5,5)+8*C(1,3)*C(5,5)+4*C(3,3)*C(5,5);  %checked many times   
    b=-4*C(1,3)^2+4*C(1,1)*C(3,3)-8*C(1,3)*C(5,5)-8*C(3,3)*C(5,5)-2*C(3,3)*invhs-2*C(5,5)*invhs;
    c=4*C(3,3)*C(5,5);
    quard1=sqrt((-b+sqrt(b^2-4*a*c))/(2*a));    
    quard2=sqrt((-b-sqrt(b^2-4*a*c))/(2*a));
    if w==1
        mu=quard1;  %mu=sin(thetap) 
    elseif w==2
        mu=quard2;
    else
        mu=sqrt(hs^2*C(5,5)/(1-hs^2*C(6,6)+hs^2*C(5,5)));
    end
end
return
