subroutine read_input(const)


    use AllType
    use constant
    
    implicit none
    

    
    type(const_struct) :: const

    character*80 :: EchoLine
    integer :: i


    
    open(10, file='TRV80_2rot.txt')
    
!    !     skip next five lines
!    do i=1,5
!        read (10, *) EchoLine
!        !write (*,*) EchoLine
!    end do

    ! aircraft weight
    read(10, 9040) const%W
    
    !write(*,*) const%W
    
    ! number of rotors [ND]
    read(10, '(1x,I12)') const%nrot
    
    !write(*,*) const%nRot

    ! number of blades for each rotor [ND]
    read(10, 9040) const%Nb(1)
    read(10, 9040) const%Nb(2)
    
    !write(*,*) const%Nb

    ! rotor radius [ft] (for now assum equal radius across all rotors)
    read(10, 9040) const%R(1)
    read(10, 9040) const%R(2)
    
    !write(*,*) const%R
    
    ! blade chord
    read(10, 9040) const%c(1)
    read(10, 9040) const%c(2)
    
    !write(*,*) const%c
    
    ! rotor angular speed [rad/s]
    read(10, 9040) const%omega(1)
    read(10, 9040) const%omega(2)
    
    !write(*,*) const%omega
    
    ! rotor verse of rotation [ND]: 1: CCW, -1: CCW
    read(10, 9040) const%sRot(1)
    read(10, 9040) const%sRot(2)
    
    !write(*,*) const%sRot
    
    ! twist angle for each rotor [rad]
    read(10, 9040) const%twist(1)
    read(10, 9040) const%twist(2)
    
    const%twist(1) = const%twist(1)*(pi/180)
    const%twist(2) = const%twist(2)*(pi/180)
    
    
    
    ! rotor solidity [ND]
    do i = 1,2
        const%sigma(i) = (const%Nb(i)*const%c(i)*pi)/const%R(i)
    end do
    
    
    ! rotor blade lift curve slope [1/rad]
    read(10, 9040) const%a0(1)
    read(10, 9040) const%a0(2)
    
    
    ! blade drag coefficients [ND]
    read(10, 9040) const%delta0(1)
    read(10, 9040) const%delta0(2)
    
    
    read(10, 9040) const%delta2(1)
    read(10, 9040) const%delta2(2)
    
    
    ! rotors position in body axes (i.e., with respect to CG) [ft]
    ! (if changing any of these, re-run G matrix computation)
    read(10, 9040) const%xRot(1)
    read(10, 9040) const%xRot(2)
    !write(*,*) const%xRot
    

    read(10, 9040) const%yRot(1)
    read(10, 9040) const%yRot(2)
    
    
    read(10, 9040) const%zRot(1)
    read(10, 9040) const%zRot(2)
    
    
    
9040 format (1x, f12.6)

    end subroutine