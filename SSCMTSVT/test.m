clear
clc
% number of radial elements [ND]
const.nRE=30; 
% number of azimuthal elements [ND]
const.nAE=360; 
% radial increments [ND]
const.dr=1/const.nRE; 
% azimutal increments [rad]
const.dpsi=2*pi/const.nAE; 

const.rSeg=[const.dr/2:const.dr:1-const.dr/2];
% vector of rotor azimuths 
const.psiv=[0:const.dpsi:2*pi-const.dpsi]; 

const.chiv=[0:10:90]*pi/180;

innerInt1 = [1,2,3;4,5,6,;7,8,9]

K=bsvel(0.9832,-0.0172,-0.1588,const.psiv,const.chiv(10));

%innerInt1 = trapz(const.psiv,K)
Gtemp = zeros(3,3,10)

Gtemp(:,:,1) = [0.3430 -0.0011 0.0000; -0.0028 0.0771 0.0000; 0.0000 0.0000 0.0779]

interp1(const.chiv,Gtemp,0)
