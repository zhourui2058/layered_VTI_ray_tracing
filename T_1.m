%Matlab two-point ray tracing forward modeling for horizontlly layered VTI media. The code is good for forward modeling with realistic earth material.
%It is not optimal for inversion as special treatments need to be applied for inversion.
%The code is a modification of the ray tracing code direct1.f from the double-difference location method at http://www.ldeo.columbia.edu/~felixw/HYPODD/HYPODD_1.3.tar.gz
%Interface processing is rigorously taken care of compare to the horizontal slowness shooting ray tracing method (https://library.seg.org/doi/abs/10.1190/geo2019-0325.1) 
%Please report any bugs and improvements to: Rui Zhou, zhourui2058@gmail.com, zhourui2058@126.com, https://www.linkedin.com/in/zhourui2058/,

function [tdir,trip]=T_1(C_all_L,top_all,rec,ms,ir,im)
loc_2d3d='2D';

nC=5;
nl_all=length(top_all)-1;
Cnl_all=nC*nl_all;

if strcmp(loc_2d3d,'2D')==1
    recr_in=rec(ir,1);
    recz_in=rec(ir,2);
    msr_in=ms(im,1);
    msz_in=ms(im,2); %2 instead of 3 
    h_dist=abs(msr_in-recr_in); %ms(im,2)=rec(ir,2)=0
else %3D
    recz_in=rec(ir,3);
    msz_in=ms(im,3);
    h_dist=sqrt((ms(im,1)-rec(ir,1))^2+(ms(im,2)-rec(ir,2))^2); 
end

pC=reshape(C_all_L(1:Cnl_all),nC,nl_all);  
C_all=zeros(6,6,nl_all);
for i=1:nl_all
    C_all(:,:,i)=[pC(1,i),pC(1,i)-2*pC(5,i),pC(2,i),0,0,0;  %pC(1,1,i),pC(1,3,i),pC(3,3,i),pC(5,5,i),pC(6,6)
                    pC(1,i)-2*pC(5,i),pC(1,i),pC(2,i),0,0,0;
                    pC(2,i),          pC(2,i),pC(3,i),0,0,0;
                    0,0,0,                      pC(4,i),0,0;
                    0,0,0,                      0,pC(4,i),0;
                    0,0,0,                      0,0,pC(5,i);];
end

[C,top,depth,~,~]=vtopdepthrecal(C_all,top_all,msz_in,recz_in);  %recalculate with RECEIVER DEPTH AT ZERO, horizontal coordinate does NOT change

if (length(top)==2 && depth==0) %single layer and ms is the same depth as rec, slight adjustment to avoid NaN
    depth=eps;
end

nl=length(top)-1;
if h_dist==0
    h_dist=eps;
end

d=zeros(1,nl);           %thickness of the layers 
for i=1:nl-1
    d(i)=top(i+1)-top(i);
end
d(nl)=depth-top(nl);     %distance of ms depth to the top of the layer that contains the ms

tdir=zeros(3,1);V=zeros(3,nl);G=V;seg=V;thg=V;th=V;n=zeros(2,3,nl);r=n;

thomsen=zeros(nC,nl);F=zeros(3,3,nl);M=zeros(7,3,nl);
for i=1:nl
    thomsen(:,i)=C2thomsen(C(:,:,i));
    [F(:,:,i), M(:,:,i)]=setMat(thomsen(:,i));
end

if nl>1  %15 lines for determining whether all layers are equal, if so, turn into a single layer
    all_equal=0; 
    for i=1:nl-1
        if sum(sum(abs(C(:,:,i)-C(:,:,nl))))==0
            all_equal=all_equal+1;
        end
    end
    if all_equal==nl-1
        nl=1;
        depth=0;
        for i=1:nl
            depth=depth+d(i);
        end
    end
end

trip=zeros(1,nl); %only for SV wave can there be trip_lication
if nl==1
    [n_temp, r_temp, ~, G_temp]=velVTI(thomsen(:,nl), F(:,:,nl), M(:,:,nl),[h_dist/sqrt(h_dist^2+depth^2),depth/sqrt(h_dist^2+depth^2)]', 'group');
    if size(n_temp,2) == 5 %check: because size(n_temp,2) is either 3 or 5, so maybe should remove temp
        n(:,1,nl)=n_temp(:,1);
        r(:,1,nl)=r_temp(:,1);
        G(1,nl)=G_temp(1);
        
        trip(1,nl)=1;
        n(:,2,nl)=NaN;
        disp('single layer trip');
        tdir(2)=NaN;
        
        n(:,3,nl)=n_temp(:,5);
        r(:,3,nl)=r_temp(:,5);
        G(3,nl)=G_temp(5);

        for w=1:3
            if w~=2
                seg(w,nl)=d(nl)./r(2,w,nl); %rsingle(2,:,nl)-->cos(thg)
                th(w,nl)=asin(n(1,w,nl));
                thg(w,nl)=asin(r(1,w,nl));
                tdir(w)=seg(w,nl)./G(w,nl);
            end
        end
    else
        r(:,:,nl)=r_temp;
        G(:,nl)=G_temp;
        seg(:,nl)=d(nl)./r(2,:,nl); %rsingle(2,:,nl)-->cos(thg)
        tdir=seg(:,nl)./G(:,nl);
    end
else %if nl > 1
    %16 lines for computing iV90min
    Vpsv90=zeros(2,nl);Vsh90=zeros(1,nl);V90=zeros(3,nl);
    for i=1:nl
        [Vpsv90(:,i)]=V_VTI_psv(C(:,:,i),[1,0]); %checked with velVTI
        [Vsh90(i)]   =V_VTI_sh (C(:,:,i),[1,0]);
        V90(:,i)=[Vpsv90(:,i);Vsh90(i)];
    end
    
    iV90=1./V90;
    iV90min=zeros(1,3);nlimin=zeros(1,3);
    for w=1:3
        iV90min(w)=1/V90(w,nl);
        nlimin(w)=nl;
        for  i=1:nl-1
            if (iV90(w,i)<iV90min(w))
                iV90min(w)=iV90(w,i);nlimin(w)=i;
            end
        end
    end
    
    mlhda=zeros(1,3);
    nb=zeros(2,3,nl);rb=nb;
    Vb=zeros(3,nl);         %phase velocity
    mlhdb=zeros(1,3);       %ms Layer Horizontal Distance
    refb=zeros(3,nl-1);     %refraction phase angle, nl>1 here
    
    nb_temp=zeros(2,3,nl);rb_temp=zeros(2,3,nl);Vb_temp=Vb;
    n_temp=zeros(2,3,nl); r_temp=zeros(2,3,nl); V_temp=Vb; G_temp=Vb;
    ubnl=zeros(3,1);
        
    adjust=0.9999;  %if adjust is too large ubnl could be too large that leads to ray tracing fails
    for w=1:3
        if nlimin(w)==nl
            ubnl(w)=1;
        else
            ubnl(w)=refx(w,iV90min(w)*adjust,1,C(:,:,nl));  
        end
        [nb_temp(:,:,nl),rb_temp(:,:,nl), Vb_temp(:,nl), ~]=velVTI(thomsen(:,nl), F(:,:,nl), M(:,:,nl), [ubnl(w),sqrt(1-ubnl(w)^2)]', 'phase');
        nb(:,w,nl)=nb_temp(:,w,nl);
        rb(:,w,nl)=rb_temp(:,w,nl);
        Vb(w,nl)=Vb_temp(w,nl);
        mlhdb(w)=d(nl)*rb(1,w,nl)/rb(2,w,nl);
        if rb(2,w,nl)==0
            mlhdb(w)=d(nl)*rb(1,w,nl)/0.0001;
        end
    end
    
    %calculate total horizontal travel distance for all the layers
    dela=mlhda; %size: 1*3, initialization
    delb=mlhdb;
    for w=1:3
        for i=1:nl-1
            refb(w,i)=refx(w,nb(1,w,nl),Vb(w,nl),C(:,:,i));
            [~, rb_temp(:,:,i), ~, ~]=velVTI(thomsen(:,i), F(:,:,i), M(:,:,i), [refb(w,i),sqrt(1-refb(w,i)^2)]', 'phase');
            rb(:,w,i)=rb_temp(:,w,i);
            delb(w)=delb(w)+d(i)*rb(1,w,i)/rb(2,w,i);
            if rb(2,w,i)==0
                delb(w)=delb(w)+d(i)*rb(1,w,i)/0.000001;
            end
        end
    end %end of for w=1:3
    
    %initialization of calculating horizontal travel distance for the ms layer
    mlhd=zeros(1,3);vec=zeros(2,3);
    cal_h_dist=zeros(1,3);
    ref=zeros(3,nl-1);
    
    for w=1:3 
        if (delb(w) < h_dist)
            disp('can NOT obtain larger delb');
            return
        end
        xtest=1000;rayloop=0;
        raythres=2*10^(-6);
        while (abs(xtest)>raythres && (w==1 || (w==2 &&trip(1)==0) || w==3))   %loop to find the zero of cal_h_dist-h_dist
            rayloop=rayloop+1;
            %input mlhda->dela, mlhdb->delb, seek mlhd for h_dist <--> input ak->f(ak), bk->f(bk), seek ck for 0,  %bisection method: mlhd(w)=(mlhdb(w)+mlhda(w))/2 is an alternative
            mlhd(w)=mlhda(w)+(h_dist-dela(w))*(mlhdb(w)-mlhda(w))/(delb(w)-dela(w)); %the method of false position 
            vec(:,w)=[mlhd(w)/sqrt(mlhd(w)^2+d(nl)^2),d(nl)/sqrt(mlhd(w)^2+d(nl)^2)]';
            [n_temp2, r_temp2, V_temp2, G_temp2]=velVTI(thomsen(:,nl), F(:,:,nl), M(:,:,nl), vec(:,w), 'group');
            if size(n_temp2,2)==5
                if w==1
                    n(:,w,nl)=n_temp2(:,w);
                    r(:,w,nl)=r_temp2(:,w);
                    V(w,nl)=V_temp2(w);
                    G(w,nl)=G_temp2(w);
                elseif w==2
                    trip(:)=1;
                    tdir(2)=NaN;
                    disp('layer nl trip');
                    xtest=10^(-100);      %ending the while loop if triplication happens
                else %w==3
                    n(:,w,nl)=n_temp2(:,5);
                    r(:,w,nl)=r_temp2(:,5);
                    V(w,nl)=V_temp2(5);
                    G(w,nl)=G_temp2(5);
                end
                
                if w==1 || w==3
                    cal_h_dist(w)=mlhd(w); %initial horizontal distance in ms layer
                    for i=1:nl-1           %refraction phase angle, following calculation is same to previous for loop
                        ref(w,i)=refx(w,n(1,w,nl),V(w,nl),C(:,:,i));  
                        [n_temp(:,:,i), r_temp(:,:,i), V_temp(:,i), G_temp(:,i)]=velVTI(thomsen(:,i), F(:,:,i), M(:,:,i), [ref(w,i),sqrt(1-ref(w,i)^2)]', 'phase');
                        n(:,w,i)=n_temp(:,w,i);r(:,w,i)=r_temp(:,w,i);V(w,i)=V_temp(w,i);G(w,i)=G_temp(w,i);
                        cal_h_dist(w)=cal_h_dist(w)+d(i)*r(1,w,i)/r(2,w,i);
                    end
                    
                    xtest=cal_h_dist(w)-h_dist;
                    if(xtest<0)             
                        mlhda(w)=mlhd(w);
                        dela(w)=cal_h_dist(w);
                    else
                        mlhdb(w)=mlhd(w);
                        delb(w)=cal_h_dist(w);
                    end
                end
            else %size(n_temp2,2)==3
                n(:,w,nl)=n_temp2(:,w);
                r(:,w,nl)=r_temp2(:,w);
                V(w,nl)=V_temp2(w);
                G(w,nl)=G_temp2(w);
                
                cal_h_dist(w)=mlhd(w); %initial calculated horizontal distance in ms layer
                for i=1:nl-1           %refraction phase angle, following calculation is same to previous for loop
                    ref(w,i)=refx(w,n(1,w,nl),V(w,nl),C(:,:,i)); %na(1,w,nl),Va(w,nl) is CORRECT
                    [n_temp(:,:,i), r_temp(:,:,i), V_temp(:,i), G_temp(:,i)]=velVTI(thomsen(:,i), F(:,:,i), M(:,:,i), [ref(w,i),sqrt(1-ref(w,i)^2)]', 'phase');
                    n(:,w,i)=n_temp(:,w,i);r(:,w,i)=r_temp(:,w,i);V(w,i)=V_temp(w,i);G(w,i)=G_temp(w,i);
                    cal_h_dist(w)=cal_h_dist(w)+d(i)*r(1,w,i)/r(2,w,i);
                end
                
                xtest=cal_h_dist(w)-h_dist;
                if(xtest<0)       %because of this, can't write into xtest=cal_h_dist(:)
                    mlhda(w)=mlhd(w);
                    dela(w)=cal_h_dist(w);
                else
                    mlhdb(w)=mlhd(w);
                    delb(w)=cal_h_dist(w);
                end
            end %end of if size(n_temp2,2) == 5 
        end     %end of while loop
        
        if w==1 || (w==2 && trip(1)==0) || w==3  % w==1 || (w==2 && trip(1)==0) || w==3
            seg(w,nl)=sqrt(mlhd(w)^2+d(nl)^2); %mlhd is from previous while loop
            th(w,nl)=asin(n(1,w,nl));
            thg(w,nl)=asin(r(1,w,nl));
            tdir(w)=seg(w,nl)/G(w,nl);         %direct-ray travel time for ms layer 
            %calculate direct-ray travel time from 1 to nl-1 layer
            for i=1:nl-1
                seg(w,i)=d(i)/r(2,w,i);
                tdir(w)=tdir(w)+seg(w,i)/G(w,i);
            end
        end
    end %end of w loop 
end %end of if nl==1

tdir=tdir+ms(im,3); %ok for tdir<0
