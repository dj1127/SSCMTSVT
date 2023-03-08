! Biot-Savart Calculation

    subroutine bsvel(x,y,z,psiv,chiv,len_psiv,K)

! -------------------------------------------------------------------------
! DESCRIPTION 
! 
! 
! INPUT 
! - x: x position of ith rotor with respect to jth rotor coordinate frame
!      [ft] 
! - y: y position of ith rotor with respect to jth rotor coordinate frame
!      [ft] 
! - z: z position of ith rotor with respect to jth rotor coordinate frame
!      [ft] 
! - psi: azimuth angle of [rad]
! - chi: skew angle [rad]
! 
! OUTPUT 
! - K: 
! 
! -------------------------------------------------------------------------
    
    integer,intent(in) :: len_psiv
    real*8, intent(in) :: x,y,z,chiv
    real*8,dimension(len_psiv),intent(in) :: psiv
    real*8,dimension(len_psiv),intent(out) :: K
    real*8,dimension(len_psiv) :: Rc
    integer :: i, j

    ! distance [ft]
    do i=1,len_psiv
        Rc(i) = SQRT(1+(x**2)+(y**2)+(z**2)-2*(x*cos(psiv(i))+y*sin(psiv(i))))
    end do
    
    ! [1/ft]
    do j=1,len_psiv
        K(j) = (1-(x*cos(psiv(j))+y*sin(psiv(j)))+Rc(j)*sin(chiv)*cos(psiv(j)))/&
        ((RC(j)+(cos(psiv(j))-x)*sin(chiv)+z*cos(chiv))*Rc(j))
        K(j) = K(j)*(-1)
    end do
    
    end subroutine