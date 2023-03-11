 ! ------------------------- SIMULATION PARAMETERS -------------------------
 
    module sim_parameter
    
    use defineAircraftProperties
    
      
    
    ! number of radial elements [ND]
    real, parameter :: nRE = 30
    
    ! number of azimuthal elemmennts [ND]
    real, parameter :: nAE = 360
    
    
    real*8 :: XSCALE(6,1), DELXLIN(6), DELCLIN(10)
    

    
    end module