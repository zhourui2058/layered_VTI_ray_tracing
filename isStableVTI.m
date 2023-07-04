%% *isStableVTI*
% Check whether input Thomsen parameters satisfy the elastic stability conditions 

%% Copyright
%
% Copyright $\copyright$ 2012 by the Society of Exploration Geophysicists
%
% For more information, go to http://software.seg.org  
%
% You must read and accept usage terms at: http://software.seg.org/disclaimer.txt before use

%% Revision history
% Original SEG version by Vladimir Grechka, Marathon Oil Company, August 2012

%% Input

% ani  - [1, 5] vector of Thomsen anisotropy parameters arranged as
%        [Vp0, Vs0, epsilon, delta, gamma]

%% Output

% flag - [scalar] equal to 1 if Thomsen parameters satisfy the elastic stability conditions and 
%                       to 0 otherwise

%%
function [flag] = isStableVTI(ani)
%% Settings
flag = 1;

%% Calculate the stiffnesses
c33 = ani(1)^2;   c11 = c33*(1 + 2*ani(3));   
c55 = ani(2)^2;   c66 = c55*(1 + 2*ani(5));

cc = (c33 - c55)*(2*ani(4)*c33 + c33 - c55);
c13 = sqrt(cc) - c55;
if cc < 0      
    fprintf('>>> Complex-valued stiffness coefficient c13 \n');
    flag = 0;   
    return;
end

%% Check the stability conditions
if ~(isreal(c11) &&  c11 > 0)
    fprintf('>>> Violation of the elastic stability condition: c11 = %g \n', c11);
    flag = 0;   
    return;
end

if ~(isreal(c33)  &&  c33 > 0)
    fprintf('>>> Violation of the elastic stability condition: c33 = %g \n', c33);
    flag = 0;   
    return;
end

if ~(isreal(c55)   &&  c55 > 0)
    fprintf('>>> Violation of the elastic stability condition: c55 = %g \n', c55);
    flag = 0;   
    return;
end

if ~(isreal(c66)  &&  c66 > 0)
    fprintf('>>> Violation of the elastic stability condition: c66 = %g \n', c66);
    flag = 0;   
    return
end

if c13^2 > c33*(c11 - c66)
    fprintf('>>> Violation of the elastic stability condition c13^2 < c33*(c11 - c66): \n');
    flag = 0;   
    return;
end

end  % of the function