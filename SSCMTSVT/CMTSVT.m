function [xdot,y] = CMTSVT(x,inp,const)

% DESCRIPTION
% State-variable CMTSVT model.  
% 
% INPUT 
% - x: state vector 
% - inp: inputs to rotor (swashplate, rotor velocities, angular rates)
% - const: data structure with rotorcraft constants
% 
% OUTPUT 
% - xdot: state derivative vector
% - y: output variables

%--------------------------------------------------------------------------

% map constants data structure to local variables 
R=const.R; 
nRot=const.nRot; 
const.sRot; 
xRot=const.xRot; 
yRot=const.yRot; 
zRot=const.zRot; 
G=const.G; 
sigma=const.sigma; 
a0=const.a0; 
sRot=const.sRot; 
twist=const.twist; 
delta0=const.delta0; 
delta2=const.delta2; 
NCTRLS=const.NCTRLS; 
rho=const.rho; 

% upack states 
lambda_s=x(1:3*nRot);
lambda_tot=x(3*nRot+1:3*nRot*2);

% unpack inputs 
inpr=reshape(inp,NCTRLS/nRot,nRot); 
swashplate=inpr(1:3,:)*pi/180;
Omega=inpr(4,:);
u=inpr(5,:);
v=inpr(6,:);
w=inpr(7,:);
p=inpr(8,:);
q=inpr(9,:);
r=inpr(10,:);

% -------------------------- ROTOR AERODYNAMICS ---------------------------

% reshape total inflow 
lambda_tot=reshape(lambda_tot,3,nRot);
% initialize variables 
F=zeros(3,nRot); 
V=zeros(1,nRot); 
chi=zeros(1,nRot); 
VT=zeros(1,nRot); 
mu=zeros(1,nRot); 
CTa=zeros(1,nRot); 
CP=zeros(1,nRot);
% aerodynamic parameters for each rotor
for i=1:nRot
    % rotor hub velocity 
    uh=u(i)+q(i)*zRot(i)-r(i)*yRot(i);
    vh=v(i)+r(i)*xRot(i)-p(i)*zRot(i);
    wh=w(i)-q(i)*xRot(i)+p(i)*yRot(i);
    % rotor hub advance ratios 
    mux=uh/Omega(i)/R(i);
    if sRot(i)==1
        % CCW rotor 
        muy=vh/Omega(i)/R(i);
    else
        % switch sign for CW rotor 
        muy=-vh/Omega(i)/R(i);
    end
    muz=wh/Omega(i)/R(i);
    % in-plane advance ratio
    mu(i)=sqrt(mux^2+muy^2);
    % wind azimuth angle
    zeta=atan2(muy,mux);
    % store commonly used trig functions of wind azimuth angle
    czeta=cos(zeta);
    szeta=sin(zeta);
    % wind to shaft transformation
    Th2w=[1  0        0
      0  czeta  -szeta
      0  szeta  czeta];
    % tansform roll and pitch rates to wind frame and non-dimensionalize
    if sRot(i)==1
        % CCW rotor
        phw=Th2w(3,:)*[r(i); q(i); p(i)]/Omega(i);
    else
        % switch sign for CW rotor 
        phw=-Th2w(3,:)*[r(i); q(i); p(i)]/Omega(i);
    end
    qhw=Th2w(2,:)*[r(i); q(i); p(i)]/Omega(i); 
    % swashplate inputs 
    theta0=swashplate(1,i); 
    % transform swashplate inputs to wind frame 
    theta1sw=Th2w(2,:)*swashplate(:,i);
    theta1cw=Th2w(3,:)*swashplate(:,i);
    % transform total inflow to wind frame 
    lambda0=lambda_tot(1,i);
    lambda1cw=Th2w(2,:)*lambda_tot(:,i);
    lambda1sw=Th2w(3,:)*lambda_tot(:,i);
    % harmonic AoA coefficients 
    alpha1sw=phw-lambda1sw+theta1sw;
    alpha1cw=qhw-lambda1cw+theta1cw;
    % thrust coefficient 
    F01=theta0*(1./3.+0.5*mu(i)^2)+0.5*mu(i)*(theta1sw+0.5*phw)+0.5*...
        (muz-lambda0)+0.25*(1.+mu(i)^2)*twist(i);
    CTa(i)=0.5*a0(i)*sigma(i)*F01;
    % set flapping angles to zero for now 
    beta0=0; 
    beta1sw=0; 
    beta1cw=0; 
    beta0_prime=0; 
    % torque and power coefficients
    delta=delta0(i)+delta2(i)*CTa(i)^2;
    CXw=0.5*a0(i)*sigma(i)*(beta0*(alpha1cw/6-0.25*mu(i)*beta0+0.125*...
        mu(i)^2*(theta1cw-beta1sw))+beta1cw*((1/6+0.125*mu(i)^2)*theta0+...
        0.125*(1+0.5*mu(i)^2)*twist(i)+0.25*(muz-lambda0)-beta0_prime/6+...
        mu(i)/16*(alpha1sw+theta1sw-3*beta1cw))+beta1sw*(mu(i)/16*...
        (alpha1cw+theta1cw-beta1sw)-0.125*mu(i)^2*beta0)-beta0_prime*...
        theta1sw/6+0.25*(muz-lambda0)*theta1sw+((alpha1sw-theta1sw)/6-...
        0.25*mu(i)*beta0_prime+0.5*mu(i)*(muz-lambda0))*theta0+((alpha1sw-...
        theta1sw)/8-mu(i)*beta0_prime/6+0.25*mu(i)*(muz-lambda0))*...
        twist(i)+0.5*(muz-lambda0)*(alpha1sw-theta1sw)-beta0_prime/3*...
        (alpha1sw-theta1sw)-0.5*mu(i)*delta/a0(i)-0.125*mu(i)^2*beta1cw*...
        (theta0+0.5*twist(i))+(0.0625*mu(i)*(alpha1cw-theta1cw)-0.125*...
        mu(i)^2*beta0-0.0625*mu(i)*beta1sw)*theta1sw+0.0625*mu(i)*...
        (alpha1sw-theta1sw)*theta1cw-0.125*mu(i)*(alpha1sw-theta1sw)*...
        beta1sw+0.375*mu(i)*(alpha1sw-theta1sw)*theta1sw);
    CP(i)=0.5*a0(i)*sigma(i)*(-(muz-lambda0)*(2*CTa(i)/(a0(i)*...
        sigma(i)))+mu(i)*(2*CXw/(a0(i)*sigma(i)))+0.25*delta/a0(i)*...
        (1.+4.6*mu(i)^2));
    % erodynamic roll and pitch moment coefficients
    F1s1=alpha1sw/3+mu(i)*(theta0+muz-lambda0+2*twist(i)/3);
    F1c1=alpha1cw/3;
    Cla=-0.5*a0(i)*sigma(i)*0.375*F1s1;
    CMa=-0.5*a0(i)*sigma(i)*0.375*F1c1;
    % aerodynamic forces and moments 
    F(:,i)=[CTa(i) CMa Cla]'; 
    % iteration to solve for quasi-steady lambda0
    lambda0_qs=lambda0;
    delta_lambda=10.;
    iter=0;
    while((abs(delta_lambda)>1e-9) && (iter<70))
        lambda0_qsn=CTa(i)/(2*sqrt(mu(i)^2+(lambda0_qs-muz)^2));
        delta_lambda=lambda0_qsn-lambda0_qs;
        lambda0_qs=lambda0_qs+0.5*delta_lambda;
        iter=iter+1;
    end
    % total advance ratio 
    VT(i)=sqrt(mu(i)^2+(muz-lambda0_qs)^2);
    % mass flow parameters
    V(i)=(mu(i)^2+(lambda0_qs-muz)*(2*lambda0_qs-muz))/VT(i);
    % wake skew angle
    chi(i)=pi/2-atan2(abs(lambda0_qs-muz),mu(i))
end

% --------------------------- INFLOW DYNAMICS -----------------------------

% initialize inflow matrix (contains self-induced and interference 
% components)
lambda=zeros(3,nRot,nRot);
% reshape self-induced inflow
lambda_s=reshape(lambda_s,3,nRot);
% fill inflow matrix with self-induced inflow 
for i=1:nRot
    lambda(:,i,i)=lambda_s(:,i);
end
% initialize L and R matrices 
L=zeros(3,3,nRot,nRot); 
Rmat=zeros(3,3,nRot,nRot); 
% populate block diagonal elements of L matrix 
for i=1:nRot
    for j=1:nRot
        if i==j
            % block diagonal elements of L matrix 
            L(:,:,i,j)=[0.5/VT(i) 0 0; 15*pi/64/VT(i)*tan(chi(i)/2) ...
                4*cos(chi(i))/V(i)/(1+cos(chi(i))) 0; 0 0 ...
                4/V(i)/(1+cos(chi(i)))];
        end 
    end 
end 
% initialize G matrix 
GG=zeros(3,3,nRot,nRot);
% populate L and R matrices
for i=1:nRot
    for j=1:nRot
        if i~=j
            % interpolate G matrix based on skew angle
            aaa = squeeze(G(:,:,i,j,:))
            Gtemp=permute(squeeze(G(:,:,i,j,:)),[3 1 2]);
            GG(:,:,i,j)=squeeze(interp1(const.chiv,Gtemp,chi(i))); 
            % off-diagonal elements of L matrix 
            Rmat(:,:,i,j)=1/VT(j)/(1-1.5*mu(j)^2)*[1 0 -3*mu(j); ...
                0 3*(1-1.5*mu(j)^2) 0; -1.5*mu(j) 0 3];
            L(:,:,i,j)=GG(:,:,i,j)*Rmat(:,:,i,j)*inv(L(:,:,j,j));
        end
    end
end
% populate inflow matrix with interference inflow (ignore transport 
% delays for now)
for i=1:nRot
    for j=1:nRot
        if i~=j 
            lambda(:,i,j)=L(:,:,i,j)*lambda(:,j,j);
        end
    end
end
% initialize new total inflow matrix 
lambda_tot_new=zeros(3,nRot);
% total inflow (self-induced + interference)
for i=1:nRot
    lambda_tot_new(:,i)=sum(lambda(:,i,:),3);
    % lambda_tot_new(:,i)=sum(lambda(:,:,i),2);
end
% initialize self-induced inflow dynamics 
lambda_dot=zeros(3,nRot);
% self-induced inflow dynamics
for i=1:nRot
    for j=1:nRot
        if i==j         
            A=const.M;
            b=F(:,i)-inv(L(:,:,i,i))*lambda(:,i,i);
            lambda_dot(:,i)=A\b;
        end
    end
end
% initialize system dynamics vector 
xdot=zeros(const.NSTATES,1); 
% inflow dynamics 
xdot(1:3*nRot)=reshape(lambda_dot,3*nRot,1); 
% add low pass filter of total inflow to dynamics 
xdot(3*nRot+1:6*nRot)=reshape((lambda_tot_new-lambda_tot)*200,3*nRot,1);

% -------------------------------- OUTPUT ---------------------------------

% thsust [lb]
Ta=CTa.*rho.*(pi*R.^2).*Omega.^2.*R.^2;
% torque [lb-ft]
Q=CP.*rho.*(pi*R.^2).*Omega(i).^3.*R(i).^2;
% power [hp]
Pow=CP.*rho.*(pi*R.^2).*Omega(i).^3.*R(i).^3/550;
% output  
y=[CTa';CP';Ta';Pow';Q']; 

return










