% DESCRIPTION 
% UH-60 helicopter (coax version) rotor parameters. 
% 
%--------------------------------------------------------------------------

% aircraft identifyier 
aircraft='H60';
% aircraft weight 
const.W=17000;
% number of rotors [ND]
const.nRot=2; 
% number of blades for each rotor [ND]
const.Nb=[4 4];
% rotor radius [ft] (for now assume equal radius across all rotors)
const.R=[26.83 26.83];
% blade chord 
const.c=[1.73 1.73];
% rotor angular speed [rad/s] 
const.Omega=[27 27]; 
% rotor verse of rotation [ND]: 1: CCW, -1: CCW 
const.sRot=[1 -1]; 
% twist angle for each rotor [rad]
const.twist=[-13 -13]*pi/180; 
% rotor solidity [ND]
const.sigma=const.Nb.*const.c./const.R/pi;
% rotor blade lift curve slope [1/rad]
const.a0=[5.73 5.73];
% blade drag coefficients [ND]
const.delta0=[0.01 0.01]; 
const.delta2=[350 350]; 
% rotors position in body axes (i.e., with respect to CG) [ft]
% (if changing any of these, re-run G matrix computation)
const.xRot=[0 0];
const.yRot=[0 0];
const.zRot=[0 -const.R(1)/2];