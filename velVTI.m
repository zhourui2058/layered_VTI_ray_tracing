%% *velVTI*
% Calculate the wavefront normals, rays, and the phase and group velocities for a given phase 
% or ray vector in VTI medium 
% Grechka, V., 2013. Ray-direction velocities in VTI mediaRay-direction velocities in VTI media. Geophysics, 78(1), pp.F1-F5.

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
% F    - [3, 3] matrix representing the Christoffel equation in the form
%        [1, p1^2, p1^4] * F * [1, p3^2, p3^4]' = 0, 
%        where p = [p1, 0, p3] is the slowness vector
% M    - [7, 3] matrix relating Psi = tan(group angle) to Theta = tan(phase angle) as
%        [Theta^0, Theta^1, ... Theta^6] * M * [Psi^0, Psi^1, Psi^2]' = 0
% vec  - [2, 1] vector specifying either phase or ray direction in the vertical plane
% flag - [string] flag equal to either 'phase' or 'group' to indicate whether variable 
%        'vec' is the phase or ray direction 

%% Output

% n    - [2, noWave] array of wavefront normals 
% r    - [2, noWave] array of rays
% V    - [1, noWave] array of phase velocities
% G    - [1, noWave] array of group velocities
%        noWave is the number of waves propagating along a fixed ray or wavefront normal direction

%% Comments
%
% * Quantities in the output arrays are arranged to correspond to the P-, SV-, and SH-waves
% * The length |noWave| of the output arrays is equal to 3 if |flag = 'phase'| and to 5 or 3 
%   if |flag = 'group'|, depending on whether the input ray crosses triplication on the SV
%   wavefront

%%
function [n, r, V, G] = velVTI(ani, F, M, vec, flag)
%% Settings
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
Psi = NaN(1,3);   Theta = NaN(1,5);   V = NaN(1,5);

tolEQ = 5.e-4;                  % tolerance for identifying real-valued and double roots in 
                                % [Theta^0, Theta^1, ... Theta^6] * M * [Psi^0, Psi^1, Psi^2]' = 0
tolISO = 1.e-8;                 % tolerance for isotropy
vec = vec/norm(vec);            % normalize variable 'vec'  
vec(abs(vec) < eps)   = eps;    % perturb variable 'vec' to avoid division by zero
vec(abs(vec) > 1/eps) = 1/eps;
Alpha = vec(1)/vec(2);

%% Various checks of the input
if nargin < 5
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Too few input arguments \n \n');
      error('>>> STOP');
end;

if ~(strcmp(flag, 'phase') == 1 || strcmp(flag, 'group') == 1)
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Incorrect value of variable flag = %s \n', flag);
    fprintf('>>> Its admissible values are ''phase'' or ''group'' \n \n');
      error('>>> STOP')
end;

if length(ani) ~= 5
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Array ''ani'' has incorrect length = %g \n', length(ani));
    fprintf('>>> Its correct length is 5 \n \n');
      error('>>> STOP')
end;

%% The P- and SV-waves 
if abs(ani(3)) < tolISO  &&  abs(ani(4)) < tolISO
    noWave = 2;                     % isotropy for the P- and SV-waves
    Psi(1:noWave) = Alpha;    Theta(1:noWave) = Alpha;    V = ani(1:noWave);

else                                % VTI
    % Compute coefficients of the polynomial
    % [Theta^0, Theta^1, ... Theta^6] * M * [Psi^0, Psi^1, Psi^2]' = 0
    pol = NaN(1,0);
    if strcmp(flag, 'group') == 1
        powers = NaN(1,3);
        for i=1:3;   powers(1,i) = Alpha^(i-1);   end;
        for i=1:7;   pol(i) = dot(M(i,:), powers, 2);   end;
    elseif strcmp(flag, 'phase') == 1
        powers = NaN(7,1);
        for i=1:7;   powers(i,1) = Alpha^(i-1);   end;
        for i=1:3;   pol(i) = dot(M(:,i), powers, 1);   end;
    end;

    % Solve equation [Theta^0, Theta^1, ... Theta^6] * M * [Psi^0, Psi^1, Psi^2]' = 0 for either
    % Theta or Psi
    allRoots  = roots(fliplr(pol));
    realRoots = real(allRoots(abs(imag(allRoots)) < tolEQ));      % real-valued roots

    % Find double roots for which the Christoffel equation is used to compute the phase velocities  
    doubleRoots = NaN(size(realRoots));
    for i=1:length(realRoots)-1
        for j=i+1:length(realRoots)
            if abs(realRoots(i) - realRoots(j)) < tolEQ*(abs(realRoots(i)) + abs(realRoots(j)))
                doubleRoots(i) = i;
                doubleRoots(j) = j;
            end;
        end;
    end;

    noWave = 0;     doubleRootIndex = 0;
    for i=1:length(realRoots)
        noWave = noWave + 1;
        if strcmp(flag, 'group') == 1
            Psi(noWave) = Alpha;          Theta(noWave) = realRoots(i);
        elseif strcmp(flag, 'phase') == 1
            Psi(noWave) = realRoots(i);   Theta(noWave) = Alpha;
        end;

        if isnan(doubleRoots(i)) == 1
            % Single root: compute the phase velocity as V = V(Theta, Psi) 
            p32 = (F(2,1)*Theta(noWave) - F(1,2)*Psi(noWave))/ ...      
                  (2*F(1,3)*Psi(noWave) - F(2,2)*Theta(noWave) + ...
                   F(2,2)*Psi(noWave)*Theta(noWave)^2 - 2*F(3,1)*Theta(noWave)^3);
            V(noWave) = 1/sqrt(p32*(1 + Theta(noWave)^2));              
        else
            % Double root: compute the phase velocity from the Christoffel equation
            f = 1 - (ani(2)/ani(1))^2;          % VTI
            tmp = atan(Theta(noWave));    s2 = sin(tmp)^2;    c2 = cos(tmp)^2;
            V(noWave) = ani(1)*sqrt(1 + ani(3)*s2 - f/2 + ...
                (-1)^doubleRootIndex*(f/2)*sqrt(1 + (4*s2/f)*(2*ani(4)*c2 - ani(3)*(c2 - s2)) + ...
                                                4*ani(3)^2*s2^2/f^2));
            doubleRootIndex = doubleRootIndex + 1;       
        end;
    end;            
end;

V = V(~isnan(V));       % remove extra NaN's
% Sort the roots to make the first one to correspond to the P-wave
[V, index] = sort(V, 'descend');

%% The SH-wave 
noWave = noWave + 1;    index(noWave) = noWave;
if strcmp(flag, 'group') == 1
    Psi(noWave) = Alpha;         Theta(noWave) = Alpha/(1 + 2*ani(5));
elseif strcmp(flag, 'phase') == 1
    Theta(noWave) = Alpha;       Psi(noWave) = Alpha*(1 + 2*ani(5));   
end;
tmp = atan(Theta(noWave));
V(noWave) = ani(2)*sqrt((1 + 2*ani(5))*sin(tmp)^2 + cos(tmp)^2);    % phase velocity

%% The rays, wavefront normals, and group velocities
n = NaN(2,noWave);   r = NaN(2,noWave);   G = NaN(1,noWave);
for i=1:length(V)
    if strcmp(flag, 'group') == 1
        tmp = atan(Theta(index(i)));
        n(:,i) = [sin(tmp); cos(tmp)];
        r(:,i) = vec;
    elseif strcmp(flag, 'phase') == 1
        tmp = atan(Psi(index(i)));
        n(:,i) = vec;
        r(:,i) = [sin(tmp); cos(tmp)];
    end;

    cosrn = dot(r(:,i), n(:,i));
    if cosrn < 0                                                    % check whether the sign 
        if strcmp(flag, 'group') == 1                               % needs to be flipped
            n(:,i) = -n(:,i);
        elseif strcmp(flag, 'phase') == 1
            r(:,i) = -r(:,i);
        end;
    end;
    G(i) = V(i)/abs(cosrn);                                         % group velocity
end;

end  % of the function
