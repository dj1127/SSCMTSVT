% AUTHOR
% Umberto Saetti
% Assistant Professor
% Department of Aerospace Engineering
% University of Maryland
% saetti@umd.edu
% 
% DATE
% 8/19/2022 
%
% DESCRIPTION
% Implementation of the state-space CMTSVT model from Ref. 1.  
% 
% REFERENCES 
% 1) Guner, F., "A Multirotor Inflow Model Based on Combined Momentum 
%    Theory and Simple Vortex Theory (CMTSVT) for Flight Simulations", 
%    VFS 78th Annual Forum, Fort Worth, TX, May 10-12, 2022. 
% 
% STATES             | INDICES     | DESCRIPTION
% _________________________________________________________________________
% lambda0_1          |  1          | uniform self-induced inflow [ND]
% lambda1c_1         |  2          | lat. self-induced inflow [ND]
% lambda1s_1         |  3          | lon. self-induced inflow [ND]
%     :              |  :          |             :
%     :              |  :          |             :
% lambda1s_n         |  nRot*3     | lon. self-induced inflow [ND]
% lambda0_tot_1      |  nRot*3+1   | uniform total inflow [ND]
% lambda1c_tot_1     |  nRot*3+2   | lat. total inflow [ND]
% lambda1s_tot_1     |  nRot*3+3   | lon. total inflow [ND]
%     :              |  :          |             :
%     :              |  :          |             :
% lambda1s_tot_n     |  nRot*6     | lon. total inflow [ND]
% 
% CONTROL INPUTS     | INDICES     | DESCRIPTION
% _________________________________________________________________________
% theta_0            |  1          | collective [rad]
% theta_1c           |  2          | lateral cyclic [rad]
% theta_1s           |  3          | longitudinal cyclic [rad] 
% u v w              |  4  5  6    | pedals [deg]
% p q r              |  7  8  9    | angula rates [rad]
%
% OUTPUTS            | INDICES     | DESCRIPTION
% _________________________________________________________________________
% CT_1               |  1          | thrust coefficient [ND]
%     :              |  :          |             :
%     :              |  :          |             :
% CT_n               |  nRot       | thrust coefficient [ND]
% CP_1               |  nRot+1     | power coefficient [ND]
%     :              |  :          |             :
%     :              |  :          |             :
% CP_n               |  2*nRot     | power coefficient [ND]
% T                  |  2*nRot+1   | thrust [lb]
%     :              |  :          |             :
%     :              |  :          |             :
% T                  |  3*nRot     | thrust [lb]
% P                  |  3*nRot+1   | power [hp]
%     :              |  :          |             :
%     :              |  :          |             :
% P                  |  4*nRot     | power [hp]
% Q                  |  4*nRot+1   | torque [lb-ft]
%     :              |  :          |             :
%     :              |  :          |             :
% Q                  |  5*nRot     | torque [lb-ft]
%
% -------------------------------------------------------------------------

clear all 
close all 
clc 

%--------------------------------- FLAGS ----------------------------------

% flag for rotor-on-rotor interaction 
% - 0: no interaction 
% - 1: enables interaction 
rotInt=1; 
% flag for computing/loading G matrix 
% - 0: load from file 
% - 1: compute 
runGmat=1; 
% flag for integration method for G matrix computation 
% 1: trapezioidal 
% 2: Gaussian quadrature 
intMethod=1; 
% flag for trim method
% 0: no trim 
% 1: trim inflow only 
% 2: trim inflow + desired CT by varying collective input (RPM fixed)
% 3: trim inflow + desired thrust by varying RPM (collective input fixed)
trimMethod=3; 

%------------------------- SIMULATION PARAMETERS --------------------------

% load aircraft parameters
% - H60: Sikorsky UH-60
% - TRV80: Malloy TRV-80 (all 4 pairs of counter-rotating coaxial rotors)
% - TRV80_2rot: Malloy TRV-80 (1 pair of counter-rotating coaxial rotors 
%               only)
% TRV80; 
% H60;
TRV80_2rot; 
% number of radial elements [ND]
const.nRE=30; 
% number of azimuthal elements [ND]
const.nAE=360; 
% radial increments [ND]
const.dr=1/const.nRE; 
% azimutal increments [rad]
const.dpsi=2*pi/const.nAE; 
% vector of non-dimensional radial stations [ND]
const.rSeg=[const.dr/2:const.dr:1-const.dr/2];
% vector of rotor azimuths 
const.psiv=[0:const.dpsi:2*pi-const.dpsi]; 
% air density [slug/ft^3]
const.rho=0.0023769;
% initialize relative distance between rotors 
xrel=zeros(const.nRot,const.nRot); 
yrel=zeros(const.nRot,const.nRot); 
zrel=zeros(const.nRot,const.nRot); 
% relative distance between the center of each rotor (normalized by rotor 
% radius)
for i=1:const.nRot
    for j=1:const.nRot
        if i~=j
            xrel(i,j)=(const.xRot(i)-const.xRot(j))/const.R(i);
            yrel(i,j)=(const.yRot(i)-const.yRot(j))/const.R(i);
            zrel(i,j)=(const.zRot(i)-const.zRot(j))/const.R(i);
        end
    end
end
% apparent mass matrix
const.M=diag([8/3/pi 16/45/pi 16/45/pi]);
% skew angle vector for G matrices calculation
const.chiv=[0:10:90]*pi/180;
   
% run only if rotor interactions are enabled
if rotInt==1
    % compute  G matrix
    if runGmat==1
        % initialize G matrix
        const.G=zeros(3,3,const.nRot,const.nRot,length(const.chiv));
        % sweep skew angle
        for iChi=1:length(const.chiv)
            for iRot=1:const.nRot
                for jRot=1:const.nRot
                    if iRot~=jRot
                        for jAero=1:length(const.psiv)
                            for iAero=1:length(const.rSeg)
                                % aerodynamic calculation point
                                xij=xrel(iRot,jRot)+const.rSeg(iAero)*...
                                    cos(const.psiv(jAero));
                                yij=yrel(iRot,jRot)+const.rSeg(iAero)*...
                                    sin(const.psiv(jAero));
                                zij=zrel(iRot,jRot);
                                % Biot-Savart calculations
                                K=bsvel(xij,yij,zij,const.psiv,...
                                    const.chiv(iChi));
                                % inner integrals
                                if intMethod==1
                                    innerInt1(iAero,jAero)=...
                                        trapz(const.psiv,K);
                                    innerInt2(iAero,jAero)=...
                                        trapz(const.psiv,K.*...
                                        cos(const.psiv));
                                    innerInt3(iAero,jAero)=...
                                        trapz(const.psiv,K.*...
                                        sin(const.psiv));
                                elseif intMethod==2
                                    spl=spline(const.psiv,K);
                                    innerInt1(iAero,jAero)=...
                                        diff(fnval(fnint(spl),[0 2*pi]));
                                    spl=spline(const.psiv,K.*...
                                        cos(const.psiv));
                                    innerInt2(iAero,jAero)=...
                                        diff(fnval(fnint(spl),[0 2*pi]));
                                    spl=spline(const.psiv,K.*...
                                        sin(const.psiv));
                                    innerInt3(iAero,jAero)=...
                                        diff(fnval(fnint(spl),[0 2*pi]));
                                end
                            end
                            % intermediate integral
                            if intMethod==1
                                interInt1(jAero)=trapz(const.rSeg,...
                                    innerInt1(:,jAero)'.*const.rSeg);
                                interInt2(jAero)=trapz(const.rSeg,...
                                    innerInt2(:,jAero)'.*const.rSeg);
                                interInt3(jAero)=trapz(const.rSeg,...
                                    innerInt3(:,jAero)'.*const.rSeg);
                                interInt4(jAero)=trapz(const.rSeg,...
                                    innerInt1(:,jAero)'.*const.rSeg.^2);
                                interInt5(jAero)=trapz(const.rSeg,...
                                    innerInt2(:,jAero)'.*const.rSeg.^2);
                                interInt6(jAero)=trapz(const.rSeg,...
                                    innerInt3(:,jAero)'.*const.rSeg.^2);
                            elseif intMethod==2
                                spl=spline(const.rSeg,...
                                    innerInt1(:,jAero)'.*const.rSeg);
                                interInt1(jAero)=...
                                    diff(fnval(fnint(spl),[0 1]));
                                spl=spline(const.rSeg,...
                                    innerInt2(:,jAero)'.*const.rSeg);
                                interInt2(jAero)=...
                                    diff(fnval(fnint(spl),[0 1]));
                                spl=spline(const.rSeg,...
                                    innerInt3(:,jAero)'.*const.rSeg);
                                interInt3(jAero)=...
                                    diff(fnval(fnint(spl),[0 1]));
                                spl=spline(const.rSeg,...
                                    innerInt1(:,jAero)'.*const.rSeg.^2);
                                interInt4(jAero)=...
                                    diff(fnval(fnint(spl),[0 1]));
                                spl=spline(const.rSeg,...
                                    innerInt2(:,jAero)'.*const.rSeg.^2);
                                interInt5(jAero)=...
                                    diff(fnval(fnint(spl),[0 1]));
                                spl=spline(const.rSeg,...
                                    innerInt3(:,jAero)'.*const.rSeg.^2);
                                interInt6(jAero)=...
                                    diff(fnval(fnint(spl),[0 1]));
                            end
                        end
                        % outer integral
                        if intMethod==1
                            G(1,1,iRot,jRot,iChi)=-1/4/pi^2*...
                                trapz(const.psiv,interInt1);
                            G(1,2,iRot,jRot,iChi)=-1/4/pi^2*...
                                trapz(const.psiv,interInt2);
                            G(1,3,iRot,jRot,iChi)=-1/4/pi^2*...
                                trapz(const.psiv,interInt3);
                            G(2,1,iRot,jRot,iChi)=-1/pi^2*...
                                trapz(const.psiv,interInt4.*...
                                cos(const.psiv));
                            G(2,2,iRot,jRot,iChi)=-1/pi^2*...
                                trapz(const.psiv,interInt5.*...
                                cos(const.psiv));
                            G(2,3,iRot,jRot,iChi)=-1/pi^2*...
                                trapz(const.psiv,interInt6.*...
                                cos(const.psiv));
                            G(3,1,iRot,jRot,iChi)=-1/pi^2*...
                                trapz(const.psiv,interInt4.*...
                                sin(const.psiv));
                            G(3,2,iRot,jRot,iChi)=-1/pi^2*...
                                trapz(const.psiv,interInt5.*...
                                sin(const.psiv));
                            G(3,3,iRot,jRot,iChi)=-1/pi^2*...
                                trapz(const.psiv,interInt6.*...
                                sin(const.psiv));
                        elseif intMethod==2
                            spl=spline(const.psiv,interInt1);
                            G(1,1,iRot,jRot,iChi)=-1/4/pi^2*...
                                diff(fnval(fnint(spl),[0 2*pi]));
                            spl=spline(const.psiv,interInt2);
                            G(1,2,iRot,jRot,iChi)=-1/4/pi^2*...
                                diff(fnval(fnint(spl),[0 2*pi]));
                            spl=spline(const.psiv,interInt3);
                            G(1,3,iRot,jRot,iChi)=-1/4/pi^2*...
                                diff(fnval(fnint(spl),[0 2*pi]));
                            spl=spline(const.psiv,interInt4.*...
                                cos(const.psiv));
                            G(2,1,iRot,jRot,iChi)=-1/pi^2*...
                                diff(fnval(fnint(spl),[0 2*pi]));
                            spl=spline(const.psiv,interInt5.*...
                                cos(const.psiv));
                            G(2,2,iRot,jRot,iChi)=-1/pi^2*...
                                diff(fnval(fnint(spl),[0 2*pi]));
                            spl=spline(const.psiv,interInt6.*...
                                cos(const.psiv));
                            G(2,3,iRot,jRot,iChi)=-1/pi^2*...
                                diff(fnval(fnint(spl),[0 2*pi]));
                            spl=spline(const.psiv,interInt4.*...
                                sin(const.psiv));
                            G(3,1,iRot,jRot,iChi)=-1/pi^2*...
                                diff(fnval(fnint(spl),[0 2*pi]));
                            spl=spline(const.psiv,interInt5.*...
                                sin(const.psiv));
                            G(3,2,iRot,jRot,iChi)=-1/pi^2*...
                                diff(fnval(fnint(spl),[0 2*pi]));
                            spl=spline(const.psiv,interInt6.*...
                                sin(const.psiv));
                            G(3,3,iRot,jRot,iChi)=-1/pi^2*...
                                diff(fnval(fnint(spl),[0 2*pi]));
                        end
                    end
                end
            end
        end
        % save G matrix to file
        str=['G_',aircraft,'.mat'];
        save(str,'G','-v7.3');
        % assign variable to data structure
        const.G=G;
    else
        % load G matrix from file
        str=['G_',aircraft,'.mat'];
        load(str);
        % assign variable to data structure
        const.G=G;
    end
else
    % set G matrix to zeros
    const.G=zeros(3,3,const.nRot,const.nRot,length(const.chiv));
end

% const.G(:,:,1,2,1)
% const.G(:,:,2,1,1)

% relative scaling of states (for each rotor)
XSCALE=ones(6,1); 
% assemble relative scale of states for all rotors 
const.XSCALE=[];
for i=1:const.nRot
    const.XSCALE=[const.XSCALE; XSCALE];
end 
% state perturbation size for trim and linearization (for each rotor)
DELXLIN=pi/180*ones(1,6)*0.1;
% assemble state perturbation size for all rotors 
const.DELXLIN=[];
for i=1:const.nRot
    const.DELXLIN=[const.DELXLIN DELXLIN];
end 
% control perturbation size for trim and linearization (for each rotor)
DELCLIN=[ones(1,3) pi/180 ones(1,3) pi/180*ones(1,3)]*0.1;
% assemble control perturbation size for all rotors 
const.DELCLIN=[];
for i=1:const.nRot
    const.DELCLIN=[const.DELCLIN DELCLIN];
end 
% relative scaling of states (for each rotor)
const.YSCALE=[1e-3*ones(1,const.nRot) 1e-4*ones(1,const.nRot) ...
    ones(1,const.nRot) ones(1,const.nRot) ones(1,const.nRot)]'; 
% number of states 
const.NSTATES=length(const.DELXLIN); 
% number of control inputs 
const.NCTRLS=length(const.DELCLIN); 
% number of outputs 
const.NOUT=5*const.nRot; 
% trim variables and trim targets 
if trimMethod==1
    % trim variables
    const.TRIMVARS=[1:const.NSTATES];
    % trim targets
    const.TRIMTARG=[1:const.NSTATES];
elseif trimMethod==2
    % trim variables
    const.TRIMVARS=[1:const.NSTATES];
    % add collective pitch to trim variables 
    for i=1:const.nRot
        const.TRIMVARS=[const.TRIMVARS const.NSTATES+1+(i-1)*...
            const.NCTRLS/const.nRot]; 
    end
    % trim targets
    const.TRIMTARG=[1:const.NSTATES];
    % add coefficient of thrust to trim targets 
    for i=1:const.nRot
        const.TRIMTARG=[const.TRIMTARG const.NSTATES+1+(i-1)]; 
    end
elseif trimMethod==3
    % trim variables
    const.TRIMVARS=[1:const.NSTATES];
    % add collective pitch to trim variables 
    for i=1:const.nRot
        const.TRIMVARS=[const.TRIMVARS const.NSTATES+4+(i-1)*...
            const.NCTRLS/const.nRot]; 
    end
    % trim targets
    const.TRIMTARG=[1:const.NSTATES];
    % add thrust to trim targets 
    for i=1:const.nRot
        const.TRIMTARG=[const.TRIMTARG const.NSTATES+5+(i-1)]; 
    end
else 
    % trim variables
    const.TRIMVARS=[1:const.NSTATES];
    % trim targets
    const.TRIMTARG=[1:const.NSTATES];
end

%---------------------------- TRIM SIMULATION -----------------------------

% swashplate inputs for each rotor 
swashplate0=[19.73 0 0];
% velocities of each rotor 
vel0=[0 0 0];
% angular rates of each rotor 
ang0=[0 0 0];
% initialize input vector 
u0=[];
% assemble input vector 
for i=1:const.nRot 
    u0=[u0; swashplate0'; const.Omega(i); vel0'; ang0'];
end 
% initial guess to inflow 
lambda0=[0.05 0 0];
lambdatot0=[0.05 0 0];
% initialize state vector 
x0=[];
% assemble state vector 
for i=1:const.nRot 
    x0=[x0; lambda0'; lambdatot0'];
end 
% trim targets
if trimMethod==1
    targ_des=zeros(const.NSTATES,1);
elseif trimMethod==2
    % desired thrust coefficient 
    for i=1:const.nRot
        CTdes(i)=(const.W/const.nRot)/(const.rho.*(pi*const.R(i).^2).*...
            const.Omega(i).^2.*const.R(i).^2);
    end
    targ_des=[zeros(const.NSTATES,1); CTdes'];
elseif trimMethod==3
    % desired thrust coefficient 
    for i=1:const.nRot
        Tdes(i)=(const.W/const.nRot);
    end
    targ_des=[zeros(const.NSTATES,1); Tdes'];
else
    targ_des=zeros(const.NSTATES,1);
end
% trim simulation 
[x0trim,u0trim,itrim]=TrimSim('CMTSVT',x0,u0,targ_des,const); 
% trim output 
[~,y0trim] = CMTSVT(x0trim,u0trim,const);

% trim induced inflow 
fprintf('\nINDUCED INFLOW\n')
fprintf('lambda0       = ')
for i=1:const.nRot
    fprintf('   %2.4f',x0trim(1+(i-1)*3))
end
fprintf('\n')
fprintf('lambda1c      = ')
for i=1:const.nRot
    fprintf('   %2.4f',x0trim(2+(i-1)*3))
end
fprintf('\n')
fprintf('lambda1s      = ')
for i=1:const.nRot
    fprintf('   %2.4f',x0trim(3+(i-1)*3))
end
fprintf('\n')
% trim total inflow 
fprintf('\nTOTAL INFLOW\n')
fprintf('lambda0       = ')
for i=1:const.nRot
    fprintf('   %2.4f',x0trim(const.nRot*3+1+(i-1)*3))
end
fprintf('\n')
fprintf('lambda1c      = ')
for i=1:const.nRot
    fprintf('   %2.4f',x0trim(const.nRot*3+2+(i-1)*3))
end
fprintf('\n')
fprintf('lambda1s      = ')
for i=1:const.nRot
    fprintf('   %2.4f',x0trim(const.nRot*3+3+(i-1)*3))
end
fprintf('\n')
% trim coefficients of thrust and power, thrust and power 
fprintf('\nPERFORMANCE METRICS\n')
fprintf('CT            = ')
for i=1:const.nRot
    fprintf('   %1.5f',y0trim(1+(i-1)))
end
fprintf('\n')
fprintf('CP            = ')
for i=1:const.nRot
    fprintf('   %1.5f',y0trim(const.nRot*1+1+(i-1)))
end
fprintf('\n')
fprintf('T [lb]        = ')
for i=1:const.nRot
    fprintf('   %5.2f',y0trim(const.nRot*2+1+(i-1)))
end
fprintf('\n')
fprintf('Q [lb-ft]     = ')
for i=1:const.nRot
    fprintf('   %5.2f',y0trim(const.nRot*4+1+(i-1)))
end
fprintf('\n')
fprintf('P [hp]        = ')
for i=1:const.nRot
    fprintf('   %5.2f',y0trim(const.nRot*3+1+(i-1)))
end
fprintf('\n')
% trim controls 
fprintf('\nROTOR CONTROLS\n')
fprintf('theta0 [deg]  = ')
for i=1:const.nRot
    fprintf('   %2.1f',u0trim(1+(i-1)*const.NCTRLS/const.nRot))
end
fprintf('\n')
fprintf('Omega [rad/s] = ')
for i=1:const.nRot
    fprintf('   %3.1f',u0trim(4+(i-1)*const.NCTRLS/const.nRot))
end
fprintf('\n')

return

%------------------------------- LINEARIZE --------------------------------

% linearize
[A,B]=LinSim('CMTSVT',x0trim,u0trim,const); 
% plot eigenvalues 
figure 
grid on 
hold on 
plot(real(eig(A)),imag(eig(A)),'*')
plot(real(eig(A(1:3,1:3))),imag(eig(A(1:3,1:3))),'ro')
plot(real(eig(A(4:6,4:6))),imag(eig(A(4:6,4:6))),'g<')
xlabel('Real','fontsize',14)
ylabel('Imag','fontsize',14)

% plot inflow or AoA distribution 

%% 









