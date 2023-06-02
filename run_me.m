%Written in Matlab2021a
clear;clc;close all;
ms=[80 151;90 155;95 161;100 165;150 161;160 165;165 171;170 175]; %microseimic location in 2D space, data reference: Li et al., 2013, GJI, Joint microseismic location and anisotropic tomography using differential arrival times and differential backazimuths
ms(:,3)=0; %origin time

C_all_L=[20.01105;7.18102335942054;16.4025;5.588496;7.15327488;22.821084;7.94939882314651;18.4041;6.4009;7.937116;18.4781646;5.60906442483534;13.198689;5.1984;7.277760;23.79951964;8.45932219941062;19.193161;7.198489;8.35024724];
C_all_L=C_all_L*10^6; 	          %stiffness tensor
 
number_rec_1=20; 
rec=zeros(20,2);
for i=1:number_rec_1; rec(i,2)=40+i*10; end   
for ir=number_rec_1+1:30    
   rec(ir,1)=110;
   rec(ir,2)=10*ir-130;           %receivers location
end
top_all=[0 95 135 195 500000000]; %formation depth

n_rec=size(rec,1);
n_ms=size(ms,1);    
parfor ir=1:n_rec 
    tdir_br_temp=zeros(n_ms,3);
    for im=1:n_ms
        [tdir_br_temp(im,:)]=T_1(C_all_L,top_all,rec,ms,ir,im);  
    end
    tdir_br(ir,:,:)=tdir_br_temp;
end

nmmnr=n_rec*n_ms;tdir=zeros(3*nmmnr,1);
for w=1:3
    for ir=1:n_rec
        for im=1:n_ms
            tdir(im+(ir-1)*n_ms+(w-1)*nmmnr)=tdir_br(ir,im,w);   
        end
    end
end

plot(tdir);grid minor;xlabel('Count');ylabel('Traveltime');
