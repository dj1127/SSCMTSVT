module AllType
    
    implicit none


    type :: const_struct

    !---------------------------- Aircraft Properties ---------------------------
        !------------------------------------------------------------------------
        !
        !
        ! DESCRIPTION 
        !
        ! TRV-80 rotor parameters for 1 set of counter-rotating coaxial rotors 
        ! only.
        !
        !
        !------------------------------------------------------------------------

        ! aircraft weight
        real*8 :: W
        
        ! number of rotors [ND]
        integer :: nRot

        ! number of blades for each rotor [ND]
        real*8 :: Nb(2)
        
        ! rotor radius [ft] (for now assum equal radius across all rotors)
        real*8 :: R(2)
        
        ! blade chord
        real*8 :: c(2)
        
        ! rotor anngular speed [rad/s]
        real*8 :: Omega(2)
        
        ! rotor verse of rotation [ND]: 1:CW, -1: CCW
        real*8 :: sRot(2)

        ! twist angle for each rotor [rad]
        real*8 :: twist(2)
        
        ! rotor solidity [ND]
        real*8 :: sigma(2)
        
        !rotor blade lift curve slope [1/rad]
        real*8 :: a0(2)
        
        ! blade drag coefficients [ND]
        real*8 :: delta0(2)

        real*8 :: delta2(2)

        ! rotors position in body axes (i.e., with respect to CG) [ft]
        real*8 :: xRot(2)

        real*8 :: yRot(2)

        real*8 :: zRot(2)
    
    !-------------------------- Simulation Parameters ---------------------------

        ! radial increments [ND]
        real*8 :: dr

        ! azimutal increments [rad]
        real*8 :: dpsi
        
        ! vector of non-dimensional radial stations [ND]
        real*8 :: rSeg(30)
        
        ! vector of rotor azimuths
        real*8 :: psiv(360)

        ! relative distance between rotors
        real*8 :: xrel(2,2), yrel(2,2), zrel(2,2)
        
        ! apparent mass matrix
        real*8 :: M(3,3)
        
        ! skew angle vector
        real*8 :: chiv(10)
        
        ! G matrix
        real*8 :: G(3,3,2,2,10)
        
        ! relative scaling of states
        real*8 :: XSCALE(12,1)

        real*8 :: YSCALE(10,1)

        ! state perturbation size for trim and linearization
        real*8 :: DELXLIN(12)

        ! control perturbation size for trim and linearization
        real*8 :: DELCLIN(20)
    
        ! number of states
        real*8 :: NSTATES

        ! number of control inputs
        real*8 :: NCTRLS

        ! number of outputs
        integer :: NOUT

        ! trim variables
        real*8,allocatable :: TRIMVARS(:)

        ! trim targets
        real*8,allocatable :: TRIMTARG(:)
        
    end type const_struct
    
end module AllType