%% *setMat*
% Compute matrices relating the slowness components and tangents of the group and phase angles
% of the P- and SV-waves propagating in VTI media
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

%% Output

% F    - [3, 3] matrix representing the Christoffel equation in the form
%        [1, p1^2, p1^4] * F * [1, p3^2, p3^4]' = 0, 
%        where p = [p1, 0, p3] is the slowness vector
% M    - [7, 3] matrix relating Psi = tan(group angle) to Theta = tan(phase angle) as
%        [Theta^0, Theta^1, ... Theta^6] * M * [Psi^0, Psi^1, Psi^2]' = 0

%%
function [F, M] = setMat(ani)
%% Settings
[thisFolderName, thisFileName, ~] = fileparts(mfilename('fullpath'));
F = zeros(3,3);
M = zeros(7,3);

%% Checks the input
if length(ani) ~= 5
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> Array ''ani'' has incorrect length = %g \n', length(ani));
    fprintf('>>> Its correct length = 5 \n \n');
      error('>>> STOP')
end;

%% Verify the stability conditions 
if isStableVTI(ani) == 0
    fprintf('\n>>> Function ''%s'' \n', fullfile(thisFolderName, [thisFileName, '.m']));
    fprintf('>>> The elastic stability conditions are violated \n');
    fprintf('    or stifness coefficient c13 is complex-valued \n');
    fprintf('    Please modify the velocities or anisotropy parameters \n \n');
      error('>>> STOP')
end;

%% Construct the output matrices

% Assign nonzero elements of matrix F representing the Christoffel equation in the form
% [1, p1^2, p1^4] * F * [1, p3^2, p3^4]' = 0, where p1 and p3 are the horizontal and vertical
% components of the slowness vector 
F(1,1) = 1;
F(1,2) = -ani(1)^2 - ani(2)^2;
F(1,3) = ani(1)^2*ani(2)^2;
F(2,1) = -((1 + 2*ani(3))*ani(1)^2) - ani(2)^2;
F(2,2) = 2*(-ani(4) + ani(3))*ani(1)^4 + 2*(1 + ani(4))*ani(1)^2*ani(2)^2;
F(3,1) = (1 + 2*ani(3))*ani(1)^2*ani(2)^2;

if nargout == 1;   return;   end;

% Construct matrix M that relates Psi = tan(group angle) to Theta = tan(phase angle) as
% [Theta^0, Theta^1, ... Theta^6] * M * [Psi^0, Psi^1, Psi^2]' = 0
M(1,3) = F(1,3)*(-F(1,2)^2 + 4*F(1,1)*F(1,3));
M(2,2) = F(2,2)*( F(1,2)^2 - 4*F(1,1)*F(1,3));
M(3,1) = F(1,3)*F(2,1)^2 + F(2,2)*(-(F(1,2)*F(2,1)) + F(1,1)*F(2,2));
M(3,3) = -2*F(1,2)*F(1,3)*F(2,1) + 4*F(1,1)*F(1,3)*F(2,2);
M(4,2) = 2*(-(F(1,1)*F(2,2)^2) + F(1,2)^2*F(3,1) + F(1,3)*(F(2,1)^2 - 4*F(1,1)*F(3,1)));
M(5,1) = -2*F(1,2)*F(2,1)*F(3,1) + 4*F(1,1)*F(2,2)*F(3,1);
M(5,3) = -(F(1,2)*F(2,1)*F(2,2)) + F(1,1)*F(2,2)^2 + F(1,2)^2*F(3,1);
M(6,2) = F(2,2)*( F(2,1)^2 - 4*F(1,1)*F(3,1));
M(7,1) = F(3,1)*(-F(2,1)^2 + 4*F(1,1)*F(3,1));

end  % of the function
