function [x0trim,u0trim,itrim]=TrimSim(aircraft,x0,u0,targ_des,const)

% DESCRIPTION 
% Trims simulation model. 
% 
% INPUT 
% - aircraft: name of aircraft dynamic model 
% - x0: guess to trim state vector 
% - u0: guess to trim control input vector 
% - targ_des: desired trim target vector 
% - const: data structure with aircraft and trim properties 
% 
% OUTPUT 
% - x0trim: trim state vector 
% - u0trim: trim control input vector 
% - itrim: flag for successful completion of trim (0: trim not achieved, 
%          1: trim achieved)
% 
%--------------------------------------------------------------------------

% unpack aircraft constants
XSCALE=const.XSCALE;
YSCALE=const.YSCALE; 
TRIMVARS=const.TRIMVARS;
TRIMTARG=const.TRIMTARG;
NSTATES=const.NSTATES;
NCTRLS=const.NCTRLS;
% initial state and control input guess
x0trim=x0;
u0trim=u0;
% initialize number of iterations and error
it=0;
err=100;
% tolerance on trim error
trim_tol=5e-4;
% maximum number of iterations
itmax=100;
% trim aircraft
fprintf('\nITERATION      TRIM ERROR\n')
while ((it<itmax)&&(err>trim_tol))
    it=it+1;
    [xdot0,y0]=feval(aircraft,x0trim,u0trim,const);
    temp=[xdot0;y0];
    targvec=temp(TRIMTARG);
    targ_err=targvec-targ_des;
    XYSCALE=[XSCALE;YSCALE];
    err=max(abs(targ_err)./XYSCALE(TRIMTARG));
    fprintf('\n%2.0f             %5.4f',it,err)
    if (err>trim_tol)
        [A,B,C,D]=LinSim(aircraft,x0trim,u0trim,const);
        Jac=[A B; C D];
        Jac=Jac(TRIMTARG,TRIMVARS);
        trimvec=[x0trim;u0trim];
        trimvec(TRIMVARS)=trimvec(TRIMVARS)-0.5*pinv(Jac)*targ_err;
        x0trim=trimvec(1:NSTATES);
        u0trim=trimvec(NSTATES+1:NSTATES+NCTRLS);
    end
end
% store info on whether trim was achieved or not 
if err>trim_tol
    fprintf('\n\nWarning: trim not acheived\n');
    itrim=0;
else
    fprintf('\n\nSuccessful trim\n');
    itrim=1;
end

return