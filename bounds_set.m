%up and low bounds for constrained optimization for microseismic location and tomography
%not used in ray tracing
lubp=0.5;
lb=ms_true(1:Cnl_all)'*(1-lubp); %see Grechka, 2013, Geophysical Prospecting paper about applying constraints
ub=ms_true(1:Cnl_all)'*(1+lubp);
for i=1:Cnl_all
    if ms_true(i)<0
        lb=ms_true(i)'*(1+lubp);
        ub=ms_true(i)'*(1-lubp);
    end
end

% for i=1:nl_all
%     lb((i-1)*5+2)=-Inf;  %need to distinguish between bounds for generating synthetic initial model parameters and bounds for constrained optimization
%     ub((i-1)*5+2)=Inf;
% end

for i=1:nl_all
    if lb((i-1)*5+1) < 0
        lb((i-1)*5+1)=0; %stability constraints
    end
    if lb((i-1)*5+3) < 0
        lb((i-1)*5+3)=0;  
    end
    if lb((i-1)*5+4) < 0
        lb((i-1)*5+4)=0;  
    end
    if lb((i-1)*5+5) < 0
        lb((i-1)*5+5)=0;  
    end
end

thickness=zeros(1,nl_all);
thickness(1)=top(2);
for i=2:nl_all-1
    thickness(i)=top(i+1)-top(i);
end
thickness(nl_all)=thickness(nl_all-1); %only for the following, not used elsewhere
for i=1:nl_all-1
    lb(Cnl_all+i)=top(i+1)-thickness(i)/20;
    ub(Cnl_all+i)=top(i+1)+thickness(i)/20;
end
if lb<0
    lb(Cnl_all+1)=0; %surface
end


for i=1:n_ms
       if strcmp(loc_2d3d,'2D')==1
           lb(Cnl_all+nl_all-1+i*3-2)=ms(i,1)-20/L;       %r 
           lb(Cnl_all+nl_all-1+i*3-1)=ms(i,2)-20/L;       %z
           if lb(Cnl_all+nl_all-1+i*3-1)<0                %below earth surface
               lb(Cnl_all+nl_all-1+i*3-1)=eps;
           end   
           lb(Cnl_all+nl_all-1+i*3)=-10000;               %t
           ub(Cnl_all+nl_all-1+i*3-2)=ms(i,1)+20/L;       %r
           ub(Cnl_all+nl_all-1+i*3-1)=ms(i,2)+20/L;       %z
           ub(Cnl_all+nl_all-1+i*3)=10000;                %t
        else %3D
           lb(Cnl_all+nl_all-1+i*4-3)=ms(i,1)-20/L;       %x
           lb(Cnl_all+nl_all-1+i*4-2)=ms(i,2)-20/L;       %y
           lb(Cnl_all+nl_all-1+i*4-1)=ms(i,3)-20/L;       %z
           if lb(Cnl_all+nl_all-1+i*4-1)<0                %below earth surface
               lb(Cnl_all+nl_all-1+i*4-1)=eps;              
           end
           lb(Cnl_all+nl_all-1+i*4)=-10000;               %t
           ub(Cnl_all+nl_all-1+i*4-3)=ms(i,1)+20/L;       %x
           ub(Cnl_all+nl_all-1+i*4-2)=ms(i,2)+20/L;       %y
           ub(Cnl_all+nl_all-1+i*4-1)=ms(i,3)+20/L;       %z
           ub(Cnl_all+nl_all-1+i*4)=10000;                %t
       end
end
lb=lb';ub=ub';
