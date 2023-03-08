function [A,B,C,D]=LinSim(aircraft,x0,u0,const)

% unpack constants 
DELXLIN=const.DELXLIN;
DELCLIN=const.DELCLIN;
NSTATES=const.NSTATES;
NCTRLS=const.NCTRLS;
NOUT=const.NOUT; 

% evaluate flight dynamics 
[xdot0,y0]=feval(aircraft,x0,u0,const);
% initialize system and control matrices 
A=zeros(NSTATES);
B=zeros(NSTATES,NCTRLS);
C=zeros(NOUT,NSTATES);
D=zeros(NOUT,NCTRLS);

% system matrix 
for k=1:NSTATES
    x_p=x0;
    x_p(k)=x_p(k)+DELXLIN(k);
    [xdot_p1,y_p1]=feval(aircraft,x_p,u0,const);
    x_p(k)=x_p(k)-2*DELXLIN(k);
    [xdot_p2,y_p2]=feval(aircraft,x_p,u0,const);
    A(:,k)=(xdot_p1-xdot_p2)/(2*DELXLIN(k));
    C(:,k)=(y_p1-y_p2)/(2*DELXLIN(k));
end

% control matrix 
for k=1:NCTRLS
    u_p=u0; 
    u_p(k)=u_p(k)+DELCLIN(k);
    [xdot_p1,y_p1]=feval(aircraft,x0,u_p,const);
    u_p(k)=u_p(k)-2*DELCLIN(k);
    [xdot_p2,y_p2]=feval(aircraft,x0,u_p,const);
    B(:,k)=(xdot_p1-xdot_p2)/(2*DELCLIN(k));
    D(:,k)=(y_p1-y_p2)/(2*DELCLIN(k));
end

return