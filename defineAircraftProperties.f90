    module defineAircraftProperties
    
    implicit none
    
    type :: const_struct
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
        real*8 :: nRot

        
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
        ! (if changing any of these, re-run G matrix computation)
        real*8 :: xRot(2)

        
        real*8 :: yRot(2)

        
        real*8 :: zRot(2)

        
        real*8 :: DELXLIN(12)
        
        real*8 :: DELCLIN(20)
        
        real*8 :: NSTATES
        
        real*8 :: NCTRLS
        
        real*8 :: NOUT
        
    end type const_struct
    
    end module