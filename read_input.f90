subroutine read_input(TRV80)


    use defineAircraftProperties
    use constant
    
    implicit none
    

    
    type(trv80_2rot) :: TRV80

    character*80 :: EchoLine
    integer :: i


    
    open(10, file='TRV80_2rot.txt')
    
!    !     skip next five lines
!    do i=1,5
!        read (10, *) EchoLine
!        !write (*,*) EchoLine
!    end do

    ! aircraft weight
    read(10, 9040) TRV80%W
    
    !write(*,*) TRV80%W
    
    ! number of rotors [ND]
    read(10, 9040) TRV80%nrot
    
    !write(*,*) TRV80%nRot

    ! number of blades for each rotor [ND]
    read(10, 9040) TRV80%Nb(1)
    read(10, 9040) TRV80%Nb(2)
    
    !write(*,*) TRV80%Nb

    ! rotor radius [ft] (for now assum equal radius across all rotors)
    read(10, 9040) TRV80%R(1)
    read(10, 9040) TRV80%R(2)
    
    !write(*,*) TRV80%R
    
    ! blade chord
    read(10, 9040) TRV80%c(1)
    read(10, 9040) TRV80%c(2)
    
    !write(*,*) TRV80%c
    
    ! rotor angular speed [rad/s]
    read(10, 9040) TRV80%omega(1)
    read(10, 9040) TRV80%omega(2)
    
    !write(*,*) TRV80%omega
    
    ! rotor verse of rotation [ND]: 1: CCW, -1: CCW
    read(10, 9040) TRV80%sRot(1)
    read(10, 9040) TRV80%sRot(2)
    
    !write(*,*) TRV80%sRot
    
    ! twist angle for each rotor [rad]
    read(10, 9040) TRV80%twist(1)
    read(10, 9040) TRV80%twist(2)
    
    TRV80%twist(1) = TRV80%twist(1)*(pi/180)
    TRV80%twist(2) = TRV80%twist(2)*(pi/180)
    
    
    
    ! rotor solidity [ND]
    do i = 1,2
        TRV80%sigma(i) = (TRV80%Nb(i)*TRV80%c(i)*pi)/TRV80%R(i)
    end do
    
    
    ! rotor blade lift curve slope [1/rad]
    read(10, 9040) TRV80%a0(1)
    read(10, 9040) TRV80%a0(2)
    
    
    ! blade drag coefficients [ND]
    read(10, 9040) TRV80%delta0(1)
    read(10, 9040) TRV80%delta0(2)
    
    
    read(10, 9040) TRV80%delta2(1)
    read(10, 9040) TRV80%delta2(2)
    
    
    ! rotors position in body axes (i.e., with respect to CG) [ft]
    ! (if changing any of these, re-run G matrix computation)
    read(10, 9040) TRV80%xRot(1)
    read(10, 9040) TRV80%xRot(2)
    !write(*,*) TRV80%xRot
    

    read(10, 9040) TRV80%yRot(1)
    read(10, 9040) TRV80%yRot(2)
    
    
    read(10, 9040) TRV80%zRot(1)
    read(10, 9040) TRV80%zRot(2)
    
    
    
9040 format (1x, f12.6)

    end subroutine