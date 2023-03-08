function [K] =  bsvel(x,y,z,psiv,chi)

% DESCRIPTION 
% 
% 
% INPUT 
% - x: x position of ith rotor with respect to jth rotor coordinate frame
%      [ft] 
% - y: y position of ith rotor with respect to jth rotor coordinate frame
%      [ft] 
% - z: z position of ith rotor with respect to jth rotor coordinate frame
%      [ft] 
% - psi: azimuth angle of [rad]
% - chi: skew angle [rad]
% 
% OUTPUT 
% - K: 
% 
% -------------------------------------------------------------------------

% distance [ft]
Rc=sqrt(1+x^2+y^2+z^2-2*(x*cos(psiv)+y*sin(psiv)));
% [1/ft]
K=-(1-(x*cos(psiv)+y*sin(psiv))+Rc.*sin(chi).*cos(psiv))./...
    (Rc+(cos(psiv)-x)*sin(chi)+z*cos(chi))./Rc;

return 