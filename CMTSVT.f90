subroutine CMTSVT(x,inp,const,xdot,y)
    
    ! DESCRIPTION
    ! State-variable CMTSVT model.  
    ! 
    ! INPUT 
    ! - x: state vector 
    ! - inp: inputs to rotor (swashplate, rotor velocities, angular rates)
    ! - const: data structure with rotorcraft constants
    ! 
    ! OUTPUT 
    ! - xdot: state derivative vector
    ! - y: output variables

    !--------------------------------------------------------------------------
    
    use defineAircraftProperties
    use constant
    use func
    
    
    implicit none
    
    type(const_struct),intent(in) :: const
    
    real*8,intent(in) :: x(12), inp(20)
    real*8,intent(out),allocatable :: xdot(:,:),y(:,:)
    
    real*8,allocatable :: y_temp(:,:)
    real*8, allocatable :: F(:,:), V(:,:), chi(:,:), VT(:,:), mu(:,:), CTa(:,:), CP(:,:)
    
    real*8,dimension(2) :: R,sRot,xRot,yRot,zRot
    
    integer :: nRot
    
    real*8 :: G(3,3,2,2,10)
    
    real*8,dimension(2) :: sigma,a0,twist,delta0,delta2
    integer :: NCTRLS
    

    real*8 :: lambda_s(6),lambda_s_reshape(3,2),lambda_tot(6), lambda_tot_reshape(3,2)
    
    real*8 :: inpr(10,2)
    integer :: arr_shape(2)
    
    real*8 :: swashplate(3,2)
    real*8,dimension(2) :: Omega,uu,vv,ww,pp,qq,rr

    real*8, dimension(3) :: Th2w_row3
    real*8 :: Th2w(3,3)
    
    real*8 :: uh, vh, wh, mux, muy, muz, zeta, czeta, szeta, phw, qhw, theta0, theta1sw, &
          theta1cw, lambda0, lambda1cw, lambda1sw, alpha1sw, alpha1cw, F01, beta0, &
          beta1sw, beta1cw, beta0_prime, delta, CXw, CT, CQ, CMy

    real*8 :: F1s1, F1c1, Cla, CMa, lambda0_qs, delta_lambda, iter,lambda0_qsn


    integer :: i,j
    
    
    ! ---------------------- INFLOW DYNAMICS PARAMETERS------------------------
    
    real*8, allocatable :: lambda(:,:,:), L(:,:,:,:), Rmat(:,:,:,:), GG(:,:,:,:)

    real*8, allocatable :: lambda_tot_new(:,:), lambda_dot(:,:)

    real*8, allocatable :: A(:,:), B(:,:)
    ! --------------------------------------------------------------------------
    
    ! map constants data structure to local variables 
    R = const%R
    nRot = const%nRot
    sRot = const%sRot
    xRot = const%xRot
    yRot = const%yRot
    zRot = const%zRot
    G = const%G
    sigma = const%sigma
    a0 = const%a0
    twist = const%twist
    delta0 = const%delta0
    delta2 = const%delta2
    NCTRLS = const%NCTRLS
    
    ! unpack states
    lambda_s = x(1:3*nRot)
    lambda_tot = x(3*nRot+1:3*nRot*2)
    
    ! unpack inputs
    arr_shape(1) = NCTRLS/nRot
    arr_shape(2) = nRot
    inpr = reshape(inp,arr_shape, order = (/ 2, 1 /))
    swashplate = inpr(1:3,:)*pi/180
    Omega(:)=inpr(4,:)
    uu(:)=inpr(5,:)
    vv(:)=inpr(6,:)
    ww(:)=inpr(7,:)
    pp(:)=inpr(8,:)
    qq(:)=inpr(9,:)
    rr(:)=inpr(10,:)

    ! --------------------------------------------------------------------------
    ! allocate array size
    allocate(xdot(const%NSTATES,1))
    allocate(F(3,nRot))
    allocate(V(1,nRot))
    allocate(chi(1,nRot))
    allocate(VT(1,nRot))
    allocate(mu(1,nRot))
    allocate(CTa(1,nRot))
    allocate(CP(1,nRot))

    allocate(A(size(const%M,Dim=1),size(const%M,Dim=2)))
    allocate(B(3,1))
    
    ! -------------------------- ROTOR AERODYNAMICS ---------------------------
    
    ! reshape total inflow
    lambda_tot_reshape = reshape(lambda_tot, (/3,nRot/), order = (/ 2, 1 /) )
    ! initialize variables
    F = 0
    V = 0
    chi = 0
    VT = 0
    mu = 0
    CTa = 0
    CP = 0

    ! aerodynamic parameters for each rotor
    do i = 1,nRot
        ! Rotor hub velocity
        uh = uu(i) + qq(i)*zRot(i) - rr(i)*yRot(i)
        vh = vv(i) + rr(i)*xRot(i) - pp(i)*zRot(i)
        wh = ww(i) - qq(i)*xRot(i) + pp(i)*yRot(i)
        
        ! Rotor hub advance ratios
        mux = uh / Omega(i) / R(i)
        if (sRot(i) == 1) then
            ! CCW rotor
            muy = vh / Omega(i) / R(i)
        else
            ! Switch sign for CW rotor
            muy = -vh / Omega(i) / R(i)
        end if
        muz = wh / Omega(i) / R(i)
        
        ! in-plane advance ratio
        mu(1,i) = sqrt(mux**2 + muy**2)
        ! wind azimuth angle
        zeta = atan2(muy, mux)
        ! store commonly used trig functions of wind azimuth angle
        czeta = dcos(zeta)
        szeta = dsin(zeta)
        ! wind to shaft transformation
        Th2w = reshape( (/1., 0., 0., 0., czeta, -szeta, 0., szeta, czeta/), (/3,3/), order = (/ 2, 1 /))

        ! tansform roll and pitch rates to wind frame and non-dimensionalize
        if (sRot(i) == 1) then
            ! CCW rotor
            phw = dot_product(Th2w(3,:), (/rr(i), qq(i), pp(i)/)) / Omega(i)
        else
            ! switch sign for CW rotor 
            phw = -dot_product(Th2w(3,:), (/rr(i), qq(i), pp(i)/)) / Omega(i)
        end if
        qhw = dot_product(Th2w(2,:), (/rr(i), qq(i), pp(i)/)) / Omega(i)
        ! swashplate inputs 
        theta0 = swashplate(1, i)
        ! transform swashplate inputs to wind frame 
        theta1sw = dot_product(Th2w(2,:), swashplate(:,i))
        theta1cw = dot_product(Th2w(3,:), swashplate(:,i))
        ! transform total inflow to wind frame 
        lambda0 = lambda_tot_reshape(1,i)
        lambda1cw = dot_product(Th2w(2,:), lambda_tot_reshape(:,i))
        lambda1sw = dot_product(Th2w(3,:), lambda_tot_reshape(:,i))
        ! harmonic AoA coefficients 
        alpha1sw = phw - lambda1sw + theta1sw
        alpha1cw = qhw - lambda1cw + theta1cw
        ! thrust coefficient 
        F01 = theta0*(1./3. + 0.5*mu(1,i)**2) + 0.5*mu(1,i)*(theta1sw + 0.5*phw) + 0.5* &
            (muz - lambda0) + 0.25*(1. + mu(1,i)**2)*twist(i)
        CTa(1,i) = 0.5*a0(i)*sigma(i)*F01
        ! set flapping angles to zero for now 
        beta0 = 0.
        beta1sw = 0.
        beta1cw = 0.
        beta0_prime = 0.
        ! torque and power coefficients
        delta = delta0(i) + delta2(i)*CTa(1,i)**2

        CXw = 0.5*a0(i)*sigma(i)*(beta0*(alpha1cw/6-0.25*mu(1,i)*beta0+0.125* &
            mu(1,i)**2*(theta1cw-beta1sw))+beta1cw*((1/6+0.125*mu(1,i)**2)*theta0+ &
            0.125*(1+0.5*mu(1,i)**2)*twist(i)+0.25*(muz-lambda0)-beta0_prime/6+ &
            mu(1,i)/16*(alpha1sw+theta1sw-3*beta1cw))+beta1sw*(mu(1,i)/16* &
            (alpha1cw+theta1cw-beta1sw)-0.125*mu(1,i)**2*beta0)-beta0_prime* &
            theta1sw/6+0.25*(muz-lambda0)*theta1sw+((alpha1sw-theta1sw)/6- &
            0.25*mu(1,i)*beta0_prime+0.5*mu(1,i)*(muz-lambda0))*theta0+((alpha1sw- &
            theta1sw)/8-mu(1,i)*beta0_prime/6+0.25*mu(1,i)*(muz-lambda0))* &
            twist(i)+0.5*(muz-lambda0)*(alpha1sw-theta1sw)-beta0_prime/3* &
            (alpha1sw-theta1sw)-0.5*mu(1,i)*delta/a0(i)-0.125*mu(1,i)**2*beta1cw* &
            (theta0+0.5*twist(i))+(0.0625*mu(1,i)*(alpha1cw-theta1cw)-0.125* &
            mu(1,i)**2*beta0-0.0625*mu(1,i)*beta1sw)*theta1sw+0.0625*mu(1,i)* &
            (alpha1sw-theta1sw)*theta1cw-0.125*mu(1,i)*(alpha1sw-theta1sw)* &
            beta1sw+0.375*mu(1,i)*(alpha1sw-theta1sw)*theta1sw)
        CP(1,i) = 0.5*a0(i)*sigma(i)*(-(muz-lambda0)*(2*CTa(1,i)/(a0(i)*sigma(i))) &
            +mu(1,i)*(2*CXw/(a0(i)*sigma(i)))+0.25*delta/a0(i)*(1.+4.6*mu(1,i)**2))

        ! aerodynamic roll and pitch moment coefficients
        F1s1=alpha1sw/3+mu(1,i)*(theta0+muz-lambda0+2*twist(i)/3)
        F1c1=alpha1cw/3
        Cla=-0.5*a0(i)*sigma(i)*0.375*F1s1
        CMa=-0.5*a0(i)*sigma(i)*0.375*F1c1
        ! aerodynamic forces and moments
        F(:,i)=(/CTa(1,i), CMa, Cla/)
        ! iteration to solve for quasi-steady lambda0
        lambda0_qs=lambda0
        delta_lambda=10.
        iter=0
        do while((abs(delta_lambda)>1e-9) .and. (iter<70))
            lambda0_qsn=CTa(1,i)/(2*sqrt(mu(1,i)**2+(lambda0_qs-muz)**2))
            delta_lambda=lambda0_qsn-lambda0_qs
            lambda0_qs=lambda0_qs+0.5*delta_lambda
            iter=iter+1
        end do
        ! total advance ratio
        VT(1,i)=sqrt(mu(1,i)**2+(muz-lambda0_qs)**2)
        ! mass flow parameters
        V(1,i)=(mu(1,i)**2+(lambda0_qs-muz)*(2*lambda0_qs-muz))/VT(1,i)
        ! wake skew angle
        chi(1,i)=pi/2-atan2(abs(lambda0_qs-muz),mu(1,i))
    end do

    ! --------------------------- INFLOW DYNAMICS -----------------------------

    ! initialize inflow matrix (contains self-induced and interference components)
    lambda = 0
    allocate(lambda(3,nRot,nRot))
    ! reshape self-induced inflow
    lambda_s_reshape = reshape(lambda_s, (/3,nRot/), order = (/ 2, 1 /))
    ! fill inflow matrix with self-induced inflow 
    do i=1,nRot
        lambda(1:3,i,i) = lambda_s_reshape(1:3,i)
    end do

    ! initialize L and R matrices 
    allocate(L(3,3,nRot,nRot))
    allocate(Rmat(3,3,nRot,nRot))
    L = 0
    Rmat = 0
    ! populate block diagonal elements of L matrix 
    do i=1,nRot
        do j=1,nRot
            if (i==j) then
                ! block diagonal elements of L matrix 
                L(:,:,i,j)= reshape((/ 0.5/VT(1,i), 0., 0., 15*pi/(64*VT(1,i))*dtan(chi(1,i)/2),&
                            4*dcos(chi(1,i))/(V(1,i)*(1.0+dcos(chi(1,i)))), 0.0, 0.0, 0.0,&
                            4/V(1,i)/(1.0+dcos(chi(1,i))) /), (/3,3/), order = (/ 2, 1 /))
            end if
        end do
    end do

    ! initialize G matrix 
    allocate(GG(3,3,nRot,nRot))
    GG = 0
    ! populate L and R matrices
    do i=1,nRot
        do j=1,nRot
            if (i/=j) then
                ! interpolate G matrix based on skew angle 
                Gtemp = transpose(G(:,:,i,j,:))
                GG(:,:,i,j) = interp1(const%chiv, Gtemp, chi(1,i))
                ! off-diagonal elements of L matrix 
                Rmat(:,:,i,j) = 1.0/VT(j)/(1.0-1.5*mu(1,j)**2) * &
                                 (/1.0, 0.0, -3.0*mu(1,j), &
                                  0.0, 3.0*(1.0-1.5*mu(1,j)**2), 0.0, &
                                 -1.5*mu(1,j), 0.0, 3.0/)
                L(:,:,i,j) = matmul(GG(:,:,i,j), matmul(Rmat(:,:,i,j), &
                                   matmul(inv(L(:,:,j,j)),L(:,:,i,j))))
            end if
        end do
    end do

    ! populate inflow matrix with interference inflow (ignore transport delays for now)
    do i=1,nRot
        do j=1,nRot
            if (i/=j) then
                lambda(:,i,j) = matmul(L(:,:,i,j), lambda(:,j,j))
            end if
        end do
    end do

    ! initialize new total inflow matrix 
    allocate(lambda_tot_new(3,nRot))
    lambda_tot_new = 0
    ! total inflow (self-induced + interference)
    do i=1,nRot
        lambda_tot_new(:,i) = sum(lambda(:,i,:), dim=3)
        ! lambda_tot_new(:,i)=sum(lambda(:,:,i),2);
    end do

    ! initialize self-induced inflow dynamics 
    allocate(lambda_dot(3,nRot))
    lambda_dot = 0
    ! self-induced inflow dynamics
    do i=1,nRot
        do j=1,nRot
            if (i==j) then
                A = const%M
                b = F(:,i) - matmul(inv(L(:,:,i,i)), lambda(:,i,i))
                lambda_dot(:,i) = matmul(inv(A), b)
            end if
        end do
    end do

    ! -------------------------------- OUTPUT ---------------------------------

    ! thrust [lb]
    Ta = CTa*rho*(pi*(R**2))*(Omega**2)*(R**2)
    ! torque [lb-ft]
    Q=CP*rho*(pi*(R**2))*(Omega(i)**3)*(R(i)**2);
    ! power [hp]
    Pow = CP*rho*(pi*(R**2))*(Omega(i)**3)*(R(i)**3)/550;
    ! output
    y_temp = (/ CTa, CP, Ta, Pow, Q /)
    y = transpose(y_temp)
    allocate(y((size(CTa,Dim=2)+size(CP,Dim=2)+size(Ta,Dim=2)+size(Pow,Dim=2)+size(Q,Dim=2)),1))
    ! CTa, CP, Ta, Pow, Q

        
    end subroutine