 ! ------------------------- SIMULATION PARAMETERS -------------------------
 
    module sim_parameter
    
    use defineAircraftProperties
    
      
    
    ! number of radial elements [ND]
    real, parameter :: nRE = 30
    
    ! number of azimuthal elemmennts [ND]
    real, parameter :: nAE = 360
    
    integer, parameter :: maxdim =100   ! Every array has size of 100 instead of allocatable
    
    
    ! radial increments [ND]
    real*8 :: dr, dpsi
    
    ! Size of the array has been hard coded
    real*8 :: rSeg(30), psiv(360)

    real*8 :: xrel(2,2), yrel(2,2), zrel(2,2)
    
    ! apparent mass matrix
    real*8 :: M(3,3)
    
    real*8 :: chiv(10)
    
    real*8 :: G(maxdim,maxdim,maxdim,maxdim,maxdim)
    
    real*8 :: XSCALE(6,1), DELXLIN(6), DELCLIN(10)
    
    real*8 :: conXSCALE(12,1), conYSCALE(10,1), conDELXLIN(12), conDELCLIN(20)
    
    real*8 :: conNSTATES, conNCTRLS, conNOUT

    
    end module