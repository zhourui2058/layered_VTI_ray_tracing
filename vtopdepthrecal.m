%recalculates C, top and msz for direct wave if the receiver depth is not zero
%input
%top: top of layers, from 0, such as top=[0 86.5 180 205.6 257];
%C: stiffness in each layer, such as v=[2.763 2.726 3.4589 6.4];
%msz: source depth
%recz: receiver depth

%output:
%topnew: new top of layers with RECEIVER DEPTH AT ZERO
%vnew: new velocity count from topnew
%msznew: new source depth
function [Cnew,topnew,msznew,rnl,msnl]=vtopdepthrecal(C,top,msz,recz)
%see bottom for explain following 16 lines
adjust=10^(-10); %may also be eps
if(ismember(msz,top)==1 && ismember(recz,top)==1)
    if(find(top==msz)==find(top==recz))              %case 1: source and receiver at the same layer boundary
        disp('source and receiver at the same layer boundary, check vtopdepthrecal');        
        msz=msz+adjust;
        recz=recz+adjust;
    end
    if(find(top==msz)<find(top==recz))               %case 2: source and receiver at different layer boundary
        disp('source and receiver at different layer boundary, check vtopdepthrecal');   
        msz=msz+adjust;
        recz=recz-adjust;
    end
    if(find(top==msz)>find(top==recz))               %case 2: source and receiver at different layer boundary
        disp('source and receiver at different layer boundary, check vtopdepthrecal'); 
        msz=msz-adjust;
        recz=recz+adjust;
    end
end
if(ismember(msz,top)==1 && ismember(recz,top)==0)    %case 3: source is at the layer boundary, and receiver is not at the layer boundary
    disp('source is at the layer boundary, and receiver is not at the layer boundary, check vtopdepthrecal'); 
    msz=msz-adjust;
end
if(ismember(msz,top)==0 && ismember(recz,top)==1)    %case 4: source is not at the layer boundary, and receiver is at the layer boundary
    disp('source is not at the layer boundary, and receiver is at the layer boundary, check vtopdepthrecal'); 
    recz,top %#ok<NOPRT>
    recz=recz+adjust;
end

if msz==recz
%     disp('msz==recz, check vtopdepthrecal');
%     f_mszeqrecz=1;
    msz=msz+adjust;  
end

for i=1:length(top)
    if (recz < top(i))
        rnl=i-1;%original receiver layer
        break;
    end
end

for i=1:length(top)
    if (msz < top(i))
        msnl=i-1;%original source layer
        break;
    end
end

if (recz<msz)
    Cnew=zeros(6,6,msnl-rnl+1);
    j=1;
    for i=rnl:msnl %copy v into vnew from rnl until msnl(including msnl)
        Cnew(:,:,j)=C(:,:,i);
        j=j+1;
    end
    
    topnew=zeros(1,msnl-rnl+2);
    %topnew(1)=0;%copy top into topnew, topnew(1) is always zero
    j=2;
    for i=rnl:msnl
        topnew(j)=top(i+1)-recz;
        j=j+1;
    end
    
    msznew=msz-recz;
end

if(recz>msz)
    topnew=zeros(1,rnl-msnl+2);
    Cnew=zeros(6,6,rnl-msnl+1);
    j=2;
    for i=rnl:-1:msnl
        topnew(j)=recz-top(i);
        Cnew(:,:,j-1)=C(:,:,i);
        j=j+1;
    end

    msznew=recz-msz;
end

topnew(end)=50000000;
