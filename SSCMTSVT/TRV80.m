% DESCRIPTION 
% TRV-80 rotor parameters for all 4 pairs of counter-rotating coaxial 
% rotors. 
% 
%--------------------------------------------------------------------------

% aircraft identifyier 
aircraft='TRV80';
% aircraft weight 
const.W=61.7;
% number of rotors [ND]
const.nRot=8; 
% number of blades for each rotor [ND]
const.Nb=[2 2 2 2 2 2 2 2];
% rotor radius [ft] (for now assume equal radius across all rotors)
const.R=[1.108333 1.108333 1.108333 1.108333 1.108333 1.108333 ...
    1.108333 1.108333];
% blade chord 
const.c=[0.182 0.182 0.182 0.182 0.182 0.182 0.182 0.182];
% rotor angular speed [rad/s] 
const.Omega=[418 418 418 418 418 418 418 418]; 
% rotor verse of rotation [ND]: 1: CCW, -1: CCW 
const.sRot=[1 -1 1 -1 1 -1 1 -1]; 
% twist angle for each rotor [rad]
const.twist=[-12.35 -12.35 -12.35 -12.35 -12.35 -12.35 -12.35 -12.35]*...
    pi/180; 
% rotor solidity [ND]
const.sigma=const.Nb.*const.c./const.R/pi;
% rotor blade lift curve slope [1/rad]
const.a0=[5.73 5.73 5.73 5.73 5.73 5.73 5.73 5.73];
% blade drag coefficients [ND]
const.delta0=[0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01]; 
const.delta2=[350 350 350 350 350 350 350 350]; 
% rotors position in body axes (i.e., with respect to CG) [ft]
% (if changing any of these, re-run G matrix computation)
const.xRot=[4.2732/2 4.2732/2 4.2732/2 4.2732/2 -4.2732/2 -4.2732/2 ...
    -4.2732/2 -4.2732/2];
const.yRot=[2.6281/2 2.6281/2 -2.6281/2 -2.6281/2  2.6281/2  2.6281/2 ...
    -2.6281/2 -2.6281/2];
const.zRot=[0.176/2 -0.176/2 0.176/2 -0.176/2 0.176/2 -0.176/2 ...
    0.176/2 -0.176/2];
