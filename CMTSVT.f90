subroutine CMTSVT(x,inp,NCTRLS,xdot,y)
    
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
    
    
    implicit none
    
    type(trv80_2rot) :: TRV80
    
    
    real*8,dimension(2) :: R,sRot,xRot,yRot,zRot
    
    integer :: nRot
    
    real*8 :: G(3,3,2,2,10)
    
    real*8,dimension(2) :: sigma,a0,twist,delta0,delta2
    real*8 :: NCTRLS
    
    real*8 :: x(12), inp(20)
    real*8 :: lambda_s(6),lambda_s_reshape(3,2),lambda_tot(6), lambda_tot_reshape(3,2)
    
    real*8 :: inpr(10,2)
    integer :: arr_shape(2)
    
    real*8 :: swashplate(3,2)
    real*8,dimension(2) :: Omega,uu,vv,ww,pp,qq,rr

    real*8, dimension(3) :: Th2w_row3
    real*8 :: Th2w(3,3)
    
    real*8 :: F(3,2)
    
    real*8, dimension(2) :: V, chi, VT, mu, CTa, CP
    
    real*8 :: uh, vh, wh, mux, muy, muz, zeta, czeta, szeta, phw, qhw, theta0, theta1sw, &
          theta1cw, lambda0, lambda1cw, lambda1sw, alpha1sw, alpha1cw, F01, beta0, &
          beta1sw, beta1cw, beta0_prime, delta, CXw, CT, CQ, CMy
    integer :: i,j
    
    
    ! ---------------------- INFLOW DYNAMICS PARAMETERS------------------------
    
    real*8 :: lambda(3,2,2)
    
    real*8 :: L(3,3,2,2), Rmat(3,3,2,2)
    
    real*8 :: GG(3,3,2,2)
    
    ! map constants data structure to local variables 
    R = TRV80%R
    nRot = TRV80%nRot
    sRot = TRV80%sRot
    xRot = TRV80%xRot
    yRot = TRV80%yRot
    zRot = TRV80%zRot
    sigma = TRV80%sigma
    a0 = TRV80%a0
    twist = TRV80%twist
    delta0 = TRV80%delta0
    delta2 = TRV80%delta2
    
    !unpack states
    lambda_s = x(1:3*nRot)
    lambda_tot = x(3*nRot+1:3*nRot*2)
    
    !unpack inputs
    arr_shape(1) = NCTRLS/nRot
    arr_shape(2) = nRot
    inpr = reshape(inp,arr_shape)
    swashplate = inpr(1:3,:)*pi/180
    Omega(:)=inpr(4,:)
    uu(:)=inpr(5,:)
    vv(:)=inpr(6,:)
    ww(:)=inpr(7,:)
    pp(:)=inpr(8,:)
    qq(:)=inpr(9,:)
    rr(:)=inpr(10,:)
    
    ! -------------------------- ROTOR AERODYNAMICS ---------------------------
    
    ! reshape total inflow
    lambda_tot_reshape = reshape(lambda_tot, (/3,nRot/) )
      ! ROTOR AERODYNAMICS
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
        mu(i) = sqrt(mux**2 + muy**2)
        ! wind azimuth angle
        zeta = atan2(muy, mux)
        ! store commonly used trig functions of wind azimuth angle
        czeta = dcos(zeta)
        szeta = dsin(zeta)
        ! wind to shaft transformation
        Th2w = reshape( (/1., 0., 0., 0., czeta, -szeta, 0., szeta, czeta/), (/3,3/))

        ! tansform roll and pitch rates to wind frame and non-dimensionalize
        if (sRot(i) == 1) then
            ! CCW rotor
            phw = dot_product(Th2w(3,:), [rr(i), qq(i), pp(i)]) / Omega(i)
        else
            ! switch sign for CW rotor 
            phw = -dot_product(Th2w(3,:), [rr(i), qq(i), pp(i)]) / Omega(i)
        end if
        qhw = dot_product(Th2w(2,:), [rr(i), qq(i), pp(i)]) / Omega(i)
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
        F01 = theta0*(1./3. + 0.5*mu(i)**2) + 0.5*mu(i)*(theta1sw + 0.5*phw) + 0.5* &
            (muz - lambda0) + 0.25*(1. + mu(i)**2)*twist(i)
        CTa(i) = 0.5*a0(i)*sigma(i)*F01
        ! set flapping angles to zero for now 
        beta0 = 0
        beta1sw = 0
        beta1cw = 0
        beta0_prime = 0
        ! torque and power coefficients
        delta = delta0(i) + delta2(i)*CTa(i)**2
        CXw = 0.5*a0(i)*sigma(i)*(beta0*(alpha1cw/6. - 0.25*mu(i)*beta0 + 0.125* &
            mu(i)**2*(theta1cw - beta1sw)) + beta1cw*((1./6. + 0.125*mu(i)**2)*theta0 + &
            0.125*(1. + 0.5*mu(i)**2)*twist(i) + 0.25*(muz - lambda0) - beta0_prime/6. + &
            mu(i)/16.*(alpha1sw + theta1sw - 3.*beta1cw)) + beta1sw*(mu(i)/16.* &
            (alpha1cw + theta1cw - beta1sw) - 0.125*mu(i)**2*beta0) - beta0_prime* &
            theta1sw/6. + 0.25*(muz - lambda0))

    end do

    ! --------------------------- INFLOW DYNAMICS -----------------------------

    ! initialize inflow matrix (contains self-induced and interference components)

    ! reshape self-induced inflow
    lambda_s_reshape = reshape(lambda_s, (/3,nRot/))
    ! fill inflow matrix with self-induced inflow 
    do i=1,nRot
        lambda(1:3,i,i) = lambda_s_reshape(1:3,i)
    end do

    ! initialize L and R matrices 
    
    ! populate block diagonal elements of L matrix 
    do i=1,nRot
        do j=1,nRot
            if (i==j) then
                ! block diagonal elements of L matrix 
                L(:,:,i,j)= reshape((/ 0.5/VT(i), 0.0, 0.0, 15*pi/(64*VT(i))*dtan(chi(i)/2),&
                            4*dcos(chi(i))/(V(i)*(1.0+dcos(chi(i)))), 0.0, 0.0, 0.0,&
                            4/V(i)/(1.0+dcos(chi(i))) /), (/3,3/))
            end if
        end do
    end do

    ! initialize G matrix 
    
    ! populate L and R matrices
    do i=1,nRot
        do j=1,nRot
            if (i/=j) then
                ! interpolate G matrix based on skew angle 
                Gtemp = transpose(G(:,:,i,j,:))
                GG(:,:,i,j) = interp1(const.chiv, Gtemp, chi(i))
                ! off-diagonal elements of L matrix 
                Rmat(:,:,i,j) = 1.0/VT(j)/(1.0-1.5*mu(j)**2) * &
                                 [1.0, 0.0, -3.0*mu(j); &
                                  0.0, 3.0*(1.0-1.5*mu(j)**2), 0.0; &
                                 -1.5*mu(j), 0.0, 3.0]
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
    real(kind=8) :: lambda_tot_new(3,nRot)
    ! total inflow (self-induced + interference)
    do i=1,nRot
        lambda_tot_new(:,i) = sum(lambda(:,i,:), dim=3)
        ! lambda_tot_new(:,i)=sum(lambda(:,:,i),2);
    end do

    ! initialize self-induced inflow dynamics 
    real(kind=8) :: lambda_dot(3,nRot)
    ! self-induced inflow dynamics
    do i=1,nRot
        do j=1,nRot
            if (i==j) then
                A = const.M
                b = F(:,i) - matmul(inv(L(:,:,i,i)), lambda(:,i,i))
                lambda_dot(:,i) = matmul(inv(A), b)
            end if
        end do
    end do

        
    

        
    end subroutine