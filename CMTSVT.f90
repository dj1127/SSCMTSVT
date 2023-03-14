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
    
    use AllType
    use constant
    use func
    
    
    implicit none
    
    type(const_struct),intent(in) :: const
    
    real*8,intent(in) :: x(12), inp(20)
    real*8,intent(out),allocatable :: xdot(:),y(:,:)
    
    !real*8,allocatable :: y_temp(:,:)
    real*8 :: y_temp(1,10)
    !real*8, allocatable :: F(:,:), V(:,:), chi(:,:), VT(:,:), mu(:,:), CTa(:,:), CP(:,:)
    real*8 :: F(3,2), V(1,2), chi(1,2), VT(1,2), mu(1,2), CTa(1,2), CP(1,2)
    real*8 :: R(2,1)
    
    real*8,dimension(2) :: sRot,xRot,yRot,zRot
    
    integer :: nRot
    
    real*8 :: G(3,3,2,2,10), Gtemp(10,3,3)
    
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
    
    !real*8, allocatable :: lambda(:,:,:), lambda_sum(:,:), L(:,:,:,:), Rmat(:,:,:,:), GG(:,:,:,:)
    real*8 :: lambda(3,2,2), lambda_sum(3,2), L(3,3,2,2), Rmat(3,3,2,2), GG(3,3,2,2)

    !real*8, allocatable :: lambda_tot_new(:,:), lambda_dot(:,:)
    real*8 :: lambda_tot_new(3,2), lambda_dot(3,2)
    
    real*8 :: lambda_tot_calc(3,2)

    !real*8, allocatable :: A(:,:), b(:,:)
    real*8 :: A(3,3), b(3,1)

    real*8 :: interp_G(1,3,3)

    ! -------------------------------- Result ----------------------------------
    
    !real*8, allocatable :: Ta(:,:), Q(:,:), Pow(:,:)
    real*8 :: Ta(1,2), Q(1,2), Pow(1,2)

    ! --------------------------------------------------------------------------
    
    ! map constants data structure to local variables 
    R(:,1) = const%R
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

    !write(*,*) R
    
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
    !allocate(xdot(const%NSTATES,1))
    !allocate(F(3,nRot)) ! 3x2
    !allocate(V(1,nRot)) ! 1x2
    !allocate(chi(1,nRot))   ! 1x2
    !allocate(VT(1,nRot))    ! 1x2
    !allocate(mu(1,nRot))    ! 1x2
    !allocate(CTa(1,nRot))   ! 1x2
    !allocate(CP(1,nRot))    ! 1x2

    !allocate(A(size(const%M,Dim=1),size(const%M,Dim=2)))    ! 3x3
    !allocate(B(3,1))
    
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
        mux = uh / Omega(i) / R(i,1)
        if (sRot(i) == 1) then
            ! CCW rotor
            muy = vh / Omega(i) / R(i,1)
        else
            ! Switch sign for CW rotor
            muy = -vh / Omega(i) / R(i,1)
        end if
        muz = wh / Omega(i) / R(i,1)
        
        ! in-plane advance ratio
        mu(1,i) = sqrt(mux**2 + muy**2)
        ! wind azimuth angle
        zeta = atan2(muy, mux)
        ! store commonly used trig functions of wind azimuth angle
        czeta = dcos(zeta)
        szeta = dsin(zeta)
        ! wind to shaft transformation
        Th2W(1,1) = 1.
        Th2W(1,2) = 0.
        Th2W(1,3) = 0.
        Th2W(2,1) = 0.
        Th2W(2,2) = czeta
        Th2W(2,3) = - szeta
        Th2W(3,1) = 0.
        Th2W(3,2) = szeta
        Th2W(3,3) = czeta

        !Th2w = reshape( (/1., 0., 0., 0., czeta, -szeta, 0., szeta, czeta/), (/3,3/), order = (/ 2, 1 /))

        write(*,*) Th2W

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
        F1s1 = alpha1sw/3+mu(1,i)*(theta0+muz-lambda0+2*twist(i)/3)
        F1c1 = alpha1cw/3
        Cla = -0.5*a0(i)*sigma(i)*0.375*F1s1
        CMa = -0.5*a0(i)*sigma(i)*0.375*F1c1
        ! aerodynamic forces and moments
        F(:,i) = (/CTa(1,i), CMa, Cla/)
        ! iteration to solve for quasi-steady lambda0
        lambda0_qs = lambda0
        delta_lambda = 10.
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


    ! ! --------------------------- INFLOW DYNAMICS -----------------------------

    ! ! initialize inflow matrix (contains self-induced and interference components)
    ! !allocate(lambda(3,nRot,nRot))
    ! lambda = 0
    
    ! !allocate(lambda_sum(3,nRot))
    ! ! reshape self-induced inflow
    ! lambda_s_reshape = reshape(lambda_s, (/3,nRot/), order = (/ 2, 1 /))
    ! ! fill inflow matrix with self-induced inflow 
    ! do i=1,nRot
    !     lambda(1:3,i,i) = lambda_s_reshape(1:3,i)
    ! end do

    ! ! initialize L and R matrices 
    ! !allocate(L(3,3,nRot,nRot))
    ! !allocate(Rmat(3,3,nRot,nRot))
    ! L = 0
    ! Rmat = 0
    ! ! populate block diagonal elements of L matrix 
    ! do i=1,nRot
    !     do j=1,nRot
    !         if (i==j) then
    !             ! block diagonal elements of L matrix 
    !             L(:,:,1,1) = 0.5/VT(1,i)
    !             L(:,:,1,2) = 0.
    !             L(:,:,1,3) = 0.
    !             L(:,:,2,1) = 15*pi/(64*VT(1,i))*dtan(chi(1,i)/2)
    !             L(:,:,2,2) = 4*dcos(chi(1,i))/(V(1,i)*(1.0+dcos(chi(1,i))))
    !             L(:,:,2,3) = 0.
    !             L(:,:,3,1) = 0.
    !             L(:,:,3,2) = 0.
    !             L(:,:,3,3) = 4/V(1,i)/(1.0+dcos(chi(1,i)))

    !             !L(:,:,i,j)= reshape((/ 0.5/VT(1,i), 0., 0., 15*pi/(64*VT(1,i))*dtan(chi(1,i)/2),&
    !             !            4*dcos(chi(1,i))/(V(1,i)*(1.0+dcos(chi(1,i)))), 0.0, 0.0, 0.0,&
    !             !            4/V(1,i)/(1.0+dcos(chi(1,i))) /), (/3,3/), order = (/ 2, 1 /))
    !         end if
    !     end do
    ! end do

    ! ! initialize G matrix 
    ! !allocate(GG(3,3,nRot,nRot))
    ! GG = 0
    ! ! populate L and R matrices
    ! do i=1,nRot
    !     do j=1,nRot
    !         if (i/=j) then
    !             ! interpolate G matrix based on skew angle 
    !             Gtemp = permute(squeeze(G(:,:,i,j,:)),(/3,1,2/))
    !             interp_G = interp1(const%chiv, Gtemp, chi(1,i))
    !             GG(:,:,i,j) = reshape((/interp_G(:,:,1),interp_G(:,:,2),interp_G(:,:,3)/),(/3,3/),order = (/2,1/))

    !             ! off-diagonal elements of L matrix 
    !             Rmat(:,:,1,1) = 1.
    !             Rmat(:,:,1,2) = 0.
    !             Rmat(:,:,1,3) = -3.0*mu(1,j)
    !             Rmat(:,:,2,1) = 0.
    !             Rmat(:,:,2,2) = 3.0*(1.0-1.5*mu(1,j)**2)
    !             Rmat(:,:,2,3) = 0.
    !             Rmat(:,:,3,1) = -1.5*mu(1,j)
    !             Rmat(:,:,3,2) = 0.
    !             Rmat(:,:,3,3) = 3.
                
    !             Rmat = Rmat * 1.0/(VT(1,j)*(1.0-1.5*mu(1,j)**2))

    !             !Rmat(:,:,i,j) = 1.0/(VT(1,j)*(1.0-1.5*mu(1,j)**2)) * &
    !             !                 (/ 1.0, 0.0, -3.0*mu(1,j), &
    !             !                  0.0, 3.0*(1.0-1.5*mu(1,j)**2), 0.0, &
    !             !                 -1.5*mu(1,j), 0.0, 3.0/)
    !             L(:,:,i,j) = matmul(GG(:,:,i,j), matmul(Rmat(:,:,i,j), &
    !                                matmul(inv(L(:,:,j,j)),L(:,:,i,j))))
    !         end if
    !     end do
    ! end do

    ! ! populate inflow matrix with interference inflow (ignore transport delays for now)
    ! do i=1,nRot
    !     do j=1,nRot
    !         if (i/=j) then
    !             lambda(:,i,j) = matmul(L(:,:,i,j), lambda(:,j,j))
    !         end if
    !     end do
    ! end do

    ! ! initialize new total inflow matrix 
    ! !allocate(lambda_tot_new(3,nRot))
    ! lambda_tot_new = 0
    ! ! total inflow (self-induced + interference)
    ! do i=1,nRot
    !     lambda_sum = lambda(:,i,:)
    !     lambda_tot_new(:,i) = sum(lambda_sum, dim=2)
    !     ! lambda_tot_new(:,i)=sum(lambda(:,:,i),2);
    ! end do

    ! ! initialize self-induced inflow dynamics 
    ! !allocate(lambda_dot(3,nRot))
    ! lambda_dot = 0
    ! !!allocate(b(3,1))
    ! ! self-induced inflow dynamics
    ! do i=1,nRot
    !     do j=1,nRot
    !         if (i==j) then
    !             A = const%M
    !             b(:,1) = F(:,i) - matmul(inv(L(:,:,i,i)), lambda(:,i,i))
    !             lambda_dot(:,i) = matmul(inv(A),b(:,1)) !A\b
    !         end if
    !     end do
    ! end do

    ! ! initialize system dynamics vector 
    ! !allocate(xdot(const%NSTATES,1))
    ! ! inflow dynamics - hardcoded
    ! xdot(1:3) = lambda_dot(:,1)
    ! xdot(4:6) = lambda_dot(:,2)
    ! write(*,*) 'pass'
    ! !xdot(1,1:3*nRot)=reshape(lambda_dot,(/3*nRot,1/),order=(/2,1/)); 
    ! ! add low pass filter of total inflow to dynamics 
    ! lambda_tot_calc = (lambda_tot_new-lambda_tot_reshape)*200
    ! xdot(7:9) = lambda_tot_calc(1:3,1)
    ! xdot(10:12) = lambda_tot_calc(1:3,2)

    ! ! -------------------------------- OUTPUT ---------------------------------

    ! ! thrust [lb]
    ! !allocate(Ta(size(CTa,dim =1),size(CTa,dim =2)))
    ! Ta = CTa*rho*(pi*(R(i,1)**2))*(Omega(i)**2)*(R(i,1)**2)
    
    ! ! torque [lb-ft]
    ! !allocate(Q(size(CP,dim =1),size(CP,dim =2)))
    ! Q=CP*rho*(pi*(R(i,1)**2))*(Omega(i)**3)*(R(i,1)**2)
    
    ! ! power [hp]
    ! !allocate(Pow(size(CP,dim =1),size(CP,dim =2)))
    ! Pow = CP*rho*(pi*(R(i,1)**2))*(Omega(i)**3)*(R(i,1)**3)/550
    
    ! ! output
    ! !allocate(y_temp(1,(size(CTa,Dim=2)+size(CP,Dim=2)+size(Ta,Dim=2)+size(Pow,Dim=2)+size(Q,Dim=2))))
    ! y_temp(1,:) = (/ CTa, CP, Ta, Pow, Q /)
    ! !allocate(y((size(CTa,Dim=2)+size(CP,Dim=2)+size(Ta,Dim=2)+size(Pow,Dim=2)+size(Q,Dim=2)),1))
    ! y = transpose(y_temp)
    
    ! CTa, CP, Ta, Pow, Q


    end subroutine