Program main

    ! AUTHOR
    ! Umberto Saetti
    ! Assistant Professor
    ! Department of Aerospace Engineering
    ! University of Maryland
    ! saetti@umd.edu
    !
    ! DATE
    ! 8/19/2022 
    !
    ! DESCRIPTION
    ! Implementation of the state-space CMTSVT model from Ref. 1.  
    !
    ! REFERENCES 
    ! 1) Guner, F., "A Multirotor Inflow Model Based on Combined Momentum 
    ! Theory and Simple Vortex Theory (CMTSVT) for Flight Simulations", 
    ! VFS 78th Annual Forum, Fort Worth, TX, May 10-12, 2022. 
    !
    ! STATES             | INDICES     | DESCRIPTION
    ! _________________________________________________________________________
    ! lambda0_1          |  1          | uniform self-induced inflow [ND]
    ! lambda1c_1         |  2          | lat. self-induced inflow [ND]
    ! lambda1s_1         |  3          | lon. self-induced inflow [ND]
    !     :              |  :          |             :
    !     :              |  :          |             :
    ! lambda1s_n         |  nRot*3     | lon. self-induced inflow [ND]
    ! lambda0_tot_1      |  nRot*3+1   | uniform total inflow [ND]
    ! lambda1c_tot_1     |  nRot*3+2   | lat. total inflow [ND]
    ! lambda1s_tot_1     |  nRot*3+3   | lon. total inflow [ND]
    !     :              |  :          |             :
    !     :              |  :          |             :
    ! lambda1s_tot_n     |  nRot*6     | lon. total inflow [ND]
    !
    ! CONTROL INPUTS     | INDICES     | DESCRIPTION
    ! _________________________________________________________________________
    ! theta_0            |  1          | collective [rad]
    ! theta_1c           |  2          | lateral cyclic [rad]
    ! theta_1s           |  3          | longitudinal cyclic [rad] 
    ! u v w              |  4  5  6    | pedals [deg]
    ! p q r              |  7  8  9    | angula rates [rad]
    !
    ! OUTPUTS            | INDICES     | DESCRIPTION
    ! _________________________________________________________________________
    ! CT_1               |  1          | thrust coefficient [ND]
    !     :              |  :          |             :
    !     :              |  :          |             :
    ! CT_n               |  nRot       | thrust coefficient [ND]
    ! CP_1               |  nRot+1     | power coefficient [ND]
    !     :              |  :          |             :
    !     :              |  :          |             :
    ! CP_n               |  2*nRot     | power coefficient [ND]
    ! T                  |  2*nRot+1   | thrust [lb]
    !     :              |  :          |             :
    !     :              |  :          |             :
    ! T                  |  3*nRot     | thrust [lb]
    ! P                  |  3*nRot+1   | power [hp]
    !     :              |  :          |             :
    !     :              |  :          |             :
    ! P                  |  4*nRot     | power [hp]
    ! Q                  |  4*nRot+1   | torque [lb-ft]
    !     :              |  :          |             :
    !     :              |  :          |             :
    ! Q                  |  5*nRot     | torque [lb-ft]
    !
    ! -------------------------------------------------------------------------
    



    use defineAircraftProperties
    use constant
    use Sim_Parameter
    use Trim_Parameter

    
    implicit none
    
    type(const_struct) :: const
    


    
    integer :: i,j
    integer :: iChi, iRot, jRot
    integer :: iAero, jAero
    integer :: x, y, z, row, col
    real*8 :: xij, yij, zij
    real, allocatable :: my_array(:)
    real*8 :: innerInt1(30,360), innerInt2(30,360), innerInt3(30,360)
    real*8 :: interInt1(360), interInt2(360), interInt3(360), &
    interInt4(360), interInt5(360),interInt6(360)
    integer :: a,b,c,d,e
    real*8 :: onemat(1,2) ! need to be edited
    
    

    
    !----- integration param -----
    real*8 :: res
    real*8 :: K(360), bs_K(360)
    real*8 :: innerInt1_t(1,30), innerInt2_t(1,30), innerInt3_t(1,30)
    real*8 :: innerInt1_1D(30), innerInt2_1D(30), innerInt3_1D(30)
    real*8 :: inner1_calc(30), inner2_calc(30), inner3_calc(30),&
    inner4_calc(30), inner5_calc(30), inner6_calc(30)
    real*8 :: inter4s_calc(360), inter5s_calc(360), inter6s_calc(360),&
    inter4c_calc(360), inter5c_calc(360), inter6c_calc(360)
    
    !------- Trim Parameter -------
    real*8 :: mat_z14(1,14)
    integer :: index_calc
    

    

    
    ! --------------------------------- FLAGS ----------------------------------
    !
    ! flag for rotor-on-rotor interaction 
    ! - 0: no interaction 
    ! - 1: enables interaction 
    
    integer :: rotInt = 1
    
    ! flag for computing/loading G matrix 
    ! - 0: load from file 
    ! - 1: compute 
    
    integer :: runGmat = 1
    
    ! flag for integration method for G matrix computation 
    ! 1: trapezioidal 
    ! 2: Gaussian quadrature 
    
    integer :: intMethod = 1
    
    ! flag for trim method
    ! 0: no trim 
    ! 1: trim inflow only 
    ! 2: trim inflow + desired CT by varying collective input (RPM fixed)
    ! 3: trim inflow + desired thrust by varying RPM (collective input fixed)
    
    integer :: trimMethod = 3
    
    ! ------------------------- SIMULATION PARAMETERS --------------------------
    !
    ! load aircraft parameters
    ! - H60: Sikorsky UH-60
    ! - TRV80: Malloy TRV-80 (all 4 pairs of counter-rotating coaxial rotors)
    ! - TRV80_2rot: Malloy TRV-80 (1 pair of counter-rotating coaxial rotors only)
    ! TRV80 
    ! H60
    ! TRV80_2rot
    

    

    
    call read_input(const)
    
    
    ! radial increments [ND]

    const%dr = 1/nRE

    ! azimutal increments [rad]
    const%dpsi = 2*pi/nAE

    ! vector of non-dimensional radial stations [ND]
    !
    ! Array = [( INIT + INCR * (i-1), i=1,floor((END-INIT)/INCR+1) )]
    ! INT = dr/2
    ! INCR = dr
    ! END = 1-dr/2
    !
    const%rSeg = [((const%dr/2) + const%dr * (i-1), i=1,ceiling(((1-const%dr/2)-(const%dr/2))/const%dr+1) )]
    


    ! vector of rotor azimuths
    !
    ! INT = 0
    ! INCR = dpsi
    ! END = 2*pi-dpsi
    !
    const%psiv = [(0 + const%dpsi * (i-1), i=1,ceiling(((2*pi-const%dpsi)-0)/const%dpsi+1) )]
    

    
    ! initialize relative distance between rotors 
    const%xrel(1:const%nRot,1:const%nRot) = 0
    const%yrel(1:const%nRot,1:const%nRot) = 0
    const%zrel(1:const%nRot,1:const%nRot) = 0
    

    
    ! relative distance between the center of each rotor (normalized by rotor radius)
    do i = 1,2  !const%nRot
        do j = 1,2  !const%nRot
            const%xrel(i,j) = (const%xRot(i)-const%xRot(j))/const%R(i)
            const%yrel(i,j) = (const%yRot(i)-const%yRot(j))/const%R(i)
            const%zrel(i,j) = (const%zRot(i)-const%zRot(j))/const%R(i)
        end do
    end do
    
    
    ! apparent mass matrix
    const%M(1:3,1:3) = 0
    
    const%M(1,1) = 8/(3*pi)
    const%M(2,2) = 16/(45*pi)
    const%M(3,3) = 16/(45*pi)

    
    
    ! skew angle vector for G matrices calculation
    !
    ! Array = [( INIT + INCR * (i-1), i=1,floor((END-INIT)/INCR+1) )]
    ! INIT = 0
    ! INCR = 10
    ! END = 90
    const%chiv = [( 10 * (i-1), i=1,10 )]*(pi/180)
    
    
    print *, '==================== Integration Start ======================'   
    
    
    ! run only if rotor interactions are enabled
    if (rotInt == 1) then
     ! compute  G matrix
        if (runGmat == 1) then
            ! initialize G matrix
            G(1:3,1:3,1:const%nRot,1:const%nRot,1:size(const%chiv)) = 0

            ! sweep skew angle
            do iChi = 1,size(const%chiv) ! 10
                do iRot = 1,const%nRot ! 2
                    do jRot = 1,const%nRot ! 2
                        if (iRot /= jRot) then 
                            do jAero = 1,size(const%psiv) ! 360
                                do iAero = 1,size(const%rSeg) ! 30
                                    
                                    ! aerodynamic calculation point
                                    xij = const%xrel(iRot,jRot)+const%rSeg(iAero)*dcos(const%psiv(jAero))
                                    yij = const%yrel(iRot,jRot)+const%rSeg(iAero)*dsin(const%psiv(jAero))
                                    zij = const%zrel(iRot,jRot)
                                    
                                    ! Biot-Savart calculations
                                    call bsvel(xij,yij,zij,const%psiv,const%chiv(iChi),size(const%psiv),K)
                                    
                                    ! inner integrals
                                    if (intMethod == 1) then
                                        call trapz(const%psiv,K,360,res)
                                        innerInt1(iAero,jAero) = res
                                        call trapz(const%psiv,K*dcos(const%psiv),360,res)
                                        innerInt2(iAero,jAero) = res
                                        call trapz(const%psiv,K*dsin(const%psiv),360,res)
                                        innerInt3(iAero,jAero) = res
                                        
                                    else if (intMethod == 2) then
!                                        innerInt1(iAero,jAero) = res
!
!                                        innerInt2(iAero,jAero) = res
!
!                                        innerInt3(iAero,jAero) = res

                                    end if
                                    
                                    
                                    
                                end do
                                
                                ! intermediate integral
                                if (intMethod == 1) then
                                res = 0.
                                innerInt1_t = transpose(innerInt1(1:30,jAero:jAero))
                                do i = 1,30
                                    innerInt1_1d(i) = innerInt1_t(i,1)
                                end do
                                do i = 1,30
                                    inner1_calc(i) = innerInt1_1D(i)*const%rSeg(i)
                                end do
                                call trapz(const%rSeg,inner1_calc,30,res)
                                interInt1(jAero) = res

                                innerInt2_t = transpose(innerInt2(1:30,jAero:jAero))
                                do i = 1,30
                                    innerInt2_1d(i) = innerInt2_t(i,1)
                                end do
                                do i = 1,30
                                    inner2_calc(i) = innerInt2_1D(i)*const%rSeg(i)
                                end do
                                call trapz(const%rSeg,inner2_calc,30,res)
                                interInt2(jAero) = res

                                innerInt3_t = transpose(innerInt3(1:30,jAero:jAero))
                                do i = 1,30
                                    innerInt3_1d(i) = innerInt3_t(i,1)
                                end do
                                do i = 1,30
                                    inner3_calc(i) = innerInt3_1D(i)*const%rSeg(i)
                                end do
                                call trapz(const%rSeg,inner3_calc,30,res)
                                interInt3(jAero) = res
                                
                                res = 0.
                                do i = 1,30
                                    inner4_calc(i) = innerInt1_1D(i)*(const%rSeg(i)**2)
                                end do
                                call trapz(const%rSeg,inner4_calc,30,res)
                                interInt4(jAero) = res
                                
                                res = 0.
                                do i = 1,30
                                    inner5_calc(i) = innerInt2_1D(i)*(const%rSeg(i)**2)
                                end do
                                call trapz(const%rSeg,inner5_calc,30,res)
                                interInt5(jAero) = res
                                
                                res = 0.
                                do i = 1,30
                                    inner6_calc(i) = innerInt3_1D(i)*(const%rSeg(i)**2)
                                end do
                                call trapz(const%rSeg,inner6_calc,30,res)
                                interInt6(jAero) = res
                                
                                
                                end if
                            end do
                        
                            ! outer integral
                            if (intMethod == 1) then
                            res = 0.
                            call trapz(const%psiv,interInt1,360,res)
                            G(1,1,iRot,jRot,iChi) = -(1/(4*(pi**2)))*res
                            
                            res = 0.
                            call trapz(const%psiv,interInt2,360,res)
                            G(1,2,iRot,jRot,iChi) = -(1/(4*(pi**2)))*res
                            
                            res = 0.
                            call trapz(const%psiv,interInt3,360,res)
                            G(1,3,iRot,jRot,iChi) = -(1/(4*(pi**2)))*res
                            
                            res = 0.
                            do i = 1,360
                                inter4c_calc(i) = interInt4(i)*dcos(const%psiv(i))
                            end do
                            call trapz(const%psiv,inter4c_calc,360,res)
                            G(2,1,iRot,jRot,iChi) = -(1/(pi**2))*res
                            
                            res = 0.
                            do i = 1,360
                                inter5c_calc(i) = interInt5(i)*dcos(const%psiv(i))
                            end do
                            call trapz(const%psiv,inter5c_calc,360,res)
                            G(2,2,iRot,jRot,iChi) = -(1/(pi**2))*res
                            
                            res = 0.
                            do i = 1,360
                                inter6c_calc(i) = interInt6(i)*dcos(const%psiv(i))
                            end do
                            call trapz(const%psiv,inter6c_calc,360,res)
                            G(2,3,iRot,jRot,iChi) = -(1/(pi**2))*res
                            
                            res = 0.
                            do i = 1,360
                                inter4s_calc(i) = interInt4(i)*sin(const%psiv(i))
                            end do
                            call trapz(const%psiv,inter4s_calc,360,res)
                            G(3,1,iRot,jRot,iChi) = -(1/(pi**2))*res
                            
                            do i = 1,360
                                inter5s_calc(i) = interInt5(i)*dsin(const%psiv(i))
                            end do
                            call trapz(const%psiv,inter5s_calc,360,res)
                            G(3,2,iRot,jRot,iChi) = -(1/(pi**2))*res
                            
                            do i = 1,360
                                inter6s_calc(i) = interInt6(i)*dsin(const%psiv(i))
                            end do
                            call trapz(const%psiv,inter6s_calc,360,res)
                            G(3,3,iRot,jRot,iChi) = -(1/(pi**2))*res
                            else if (intMethod == 2) then
                            
                            end if
                        end if
                    end do
                end do
            end do
            ! save G matrix to file
            open(1, file = 'G_TRV80_2rot.dat', status='new')
            do a =1,10
                do b =1,2
                    do c = 1,2
                        do d = 1,3

                                write(1,*) (G(e,d,c,b,a),e=1,3)

                        end do
                    end do
                end do
            end do
            close(1)
            
            print *, 'Successful Output'
        
        else

        end if
    else
    ! set G matrix to zeros
        G(1:3,1:3,1:const%nRot,1:const%nRot,1:size(const%chiv)) = 0
    end if
!
!
!    
!    ! relative scaling of states (for each rotor)
!    do i = 1,6
!        XSCALE(i,1) = 1
!    end do
!    ! assemble relative scale of states for all rotors - HARDCODED
!    const%XSCALE(1:6,1) = XSCALE(1:6,1)
!    const%XSCALE(7:12,1) = XSCALE(1:6,1)
!    
!    ! state perturbation size for trim and linearization (for each rotor)
!    do j = 1,6
!        DELXLIN(j) = (pi/180.)*0.1
!    end do
!    
!    ! assemble relative scale of states for all rotors 
!    const%DELXLIN(1:6) = DELXLIN(1:6)
!    const%DELXLIN(7:12) = DELXLIN(1:6)
!    
!    ! control perturbation size for trim and linearization (for each rotor)
!    DELCLIN(1:3) = 1
!    DELCLIN(4) = pi/180
!    DELCLIN(5:7) = 1
!    DELCLIN(8:10) = pi/180
!    DELCLIN = DELCLIN*0.1
!    
!    ! assemble control perturbation size for all rotors 
!    const%DELCLIN(1:10) = DELCLIN(1:10)
!    const%DELCLIN(11:20) = DELCLIN(1:10)
!    
!    ! relative scaling of states (for each rotor)
!    call ones(1,const%nRot,onemat)
!    do i = 1,2
!        const%YSCALE(i,1) = onemat(1,i)*10.**(-3)
!    end do
!    do i = 3,4
!        const%YSCALE(i,1) = onemat(1,i-2)*10.**(-4)
!    end do
!    do i = 5,6
!        const%YSCALE(i,1) = onemat(1,i-4)
!    end do
!    do i = 7,8
!        const%YSCALE(i,1) = onemat(1,i-6)
!    end do
!    do i = 9,10
!        const%YSCALE(i,1) = onemat(1,i-8)
!    end do
!    
!    ! number of states 
!    const%NSTATES = size(conDELXLIN)
!    ! number of control inputs 
!    const%NCTRLS = size(conDELCLIN)
!    ! number of outputs 
!    const%NOUT = 5*const%nRot
!    ! trim variables and trim targets 
!!    if (trimMethod ==1) then
!!    else if (trimMethod ==2) then
!!        do i = 1,const%nRot
!!        
!!        end do
!!        
!!        do i = 1,const%nRot
!!        
!!        end do
!!    else if (trimMethod ==3) then
!!        do i = 1,const%nRot
!!        
!!        end do
!!        
!!        do i = 1,const%nRot
!!        
!!        end do
!!    else
!    
!!    end if
!    
!    
!
!!   ---------------------------- TRIM SIMULATION -----------------------------
!    
!    ! swashplate inputs for each rotor
!    swashplate0 = (/ 19.73,0.,0. /)
!    ! velocities of each rotor
!    vel0 = (/0,0,0/)
!    ! angular rates of each rotor
!    ang0 = (/0,0,0/)
!    ! initialize input vector
!    u0 = 0
!    ! assemble input vector
!    u0(1:3,1) = swashplate0
!    u0(4:5,1) = const%omega(1)
!    u0(6:8,1) = vel0
!    u0(9:11,1) = ang0
!    
!    do i = 1,3
!        u0(i,1) = swashplate0(i)
!    end do
!    do i = 4,5
!        u0(i,1) = const%omega(i-3)
!    end do
!    do i = 6,8
!        u0(i,1) = vel0(i-5)
!    end do
!    do i = 9,11
!        u0(i,1) = ang0(i-8)
!    end do
!    ! initial guess to inflow
!    lambda0 = (/0.05,0.,0./)
!    lambdatot0 = (/0.05,0.,0./)
!    ! initialize state vector
!    x0 = 0
!    ! assemble state vector
!    
!    ! trim targets
!    if (trimMethod == 1) then
!        do i =1,const%NSTATES
!            targ_des(i,1) = 0
!        end do
!    elseif (trimMethod == 2) then
!        ! desired thrust coefficient
!        do i = 1,const%nRot
!            CTdes(i) = (const%W/const%nRot)/(rho*(pi*const%R(i)**2)&
!            *const%Omega(I)**2*const%R(i)**2)
!        end do
!        do i = 1,12
!            targ_des(i,1) = 0
!        end do
!        do i = 13,14
!            targ_des(i,1) = CTdes(i-12)
!        end do
!    elseif (trimMethod == 2) then
!        ! desired thrust coefficient
!        do i = 1,const%nRot
!            Tdes(i) = const%W/const%nRot
!        end do 
!        do i = 1,12
!            targ_des(i,1) = 0
!        end do
!        do i = 13,14
!            targ_des(i,1) = Tdes(i-12)
!        end do
!    else
!        call zero(const%NSTATES,1,mat_z14)
!        targ_des(1:14,1) = mat_z14(1,1:14)
!    end if
!
!
!
!    !trim indueced inflow
!    print *, 'INDUCED INFLOW'
!    do i=1,const%nRot
!        index_calc = 1+(i-1)*3
!        x0trim_val = x0trim(index_calc)
!    end do
!    print 24,'lambda0      = ', x0trim_val
!    do i=1,const%nRot
!        index_calc = 2+(i-1)*3
!        x0trim_val = x0trim(index_calc)
!    end do
!    print 24,'lambda1c      = ', x0trim_val
!    do i=1,const%nRot
!        index_calc = 3+(i-1)*3
!        x0trim_val = x0trim(index_calc)
!    end do
!    print 24,'lambda1s      = ', x0trim_val
!    print *,' '
!    
!    ! trim total inflow 
!    print *,'TOTAL INFLOW'
!    do i=1,const%nRot
!        index_calc = const%nRot*3+1+(i-1)*3
!        x0trim_val = x0trim(index_calc)
!    end do
!    print 24,'lambda0       = ', x0trim_val
!    do i=1,const%nRot
!        index_calc = const%nRot*3+2+(i-1)*3
!        x0trim_val = x0trim(index_calc)
!    end do
!    print 24,'lambda1c      = ', x0trim_val
!    print *,' '
!    print *,'lambda1s      = '
!    do i=1,const%nRot
!        index_calc = const%nRot*3+3+(i-1)*3
!        x0trim_val = x0trim(index_calc)
!    end do
!    print 24, x0trim_val
!    print *,' '
!    
!    ! trim coefficients of thrust and power, thrust and power 
!    print *,'PERFORMANCE METRICS'
!    do i=1,const%nRot
!        index_calc = 1+(i-1)
!        y0trim_val = y0trim(index_calc)
!    end do
!    print 15,'CT            = ', y0trim_val
!    do i=1,const%nRot
!        index_calc = const%nRot*1+1+(i-1)
!        y0trim_val = y0trim(index_calc)
!    end do
!    print 15,'CP            = ', y0trim_val
!    do i=1,const%nRot
!        index_calc = const%nRot*2+1+(i-1)
!        y0trim_val = y0trim(index_calc)
!    end do
!    print 52,'T [lb]        = ', y0trim_val
!    do i=1,const%nRot
!        index_calc = const%nRot*4+1+(i-1)
!        y0trim_val = y0trim(index_calc)
!    end do
!    print 52,'Q [lb-ft]     = ', y0trim_val
!    do i=1,const%nRot
!        index_calc = const%nRot*3+1+(i-1)
!        y0trim_val = y0trim(index_calc)
!    end do
!    print 52,'P [hp]        = ', y0trim_val
!    print *,' '
!    
!    ! trim controls 
!    print *,'ROTOR CONTROLS'
!    
!    do i=1,const%nRot
!        index_calc = 1+(i-1)*const%NCTRLS/const%nRot
!        u0trim_val = u0trim(index_calc)
!    end do
!    print 21,'theta0 [deg]  = ', u0trim_val
!    
!    do i=1,const%nRot
!        index_calc = 4+(i-1)*const%NCTRLS/const%nRot
!        u0trim_val = u0trim(index_calc)
!    end do
!    print 31,'Omega [rad/s] = ', u0trim_val
!    print *,' '
!        
!
!24  format(1x, f6.4)
!
!15  format(1x, f6.5)
!
!52  format(1x, f7.2)
!
!21  format(1x, f3.1)
!
!31  format(1x, f4.1)

end program main

