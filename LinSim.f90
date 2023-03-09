subroutine LinSim(x0,u0,const,A,B,C,D)

    use defineAircraftProperties
    
    implicit none
    
    type(const_),intent(in) :: const
    
    ! const
    real*8 :: DELXLIN(12), DELCLIN(20)
    integer :: NSTATES, NCTRLS, NOUT
    
    real*8, intent(in) :: x0(12,1), u0(20,1)
    real*8, intent(out) :: A(12),B(12,20),C(10,12),D(10,20)
    integer:: k
    real*8 :: x_p(12,1), xdot0(12,1), y0(10,1)
    real*8 :: xdot_p1(12,1), y_p1(10,1), xdot_p2(12,1), y_p2(10,1)
    
    real*8 :: u_p(20,1)
    
    ! unpack constants
    DELXLIN = const%DELXLIN
    DELCLIN = const%DELCLIN
    NSTATES = const%NSTATES
    NCTRLS = const%NCTRLS
    NOUT = const%NOUT
    
    ! evaluate flight dynamics
    call CMTSVT(x0,u0,const,xdot0,y0)
    ! initializei systme and control matrices
    A = 0
    B = 0
    C = 0
    D = 0
    
    ! system matrix
    do k = 1,NSTATES
        x_p = x0
        x_p(k,1) = x_p(k,1)+DELXLIN(k)
        call CMTSVT(x_p,u0,const,xdot_p1,y_p1)

        x_p(k,1)=x_p(k,1)-2*DELXLIN(k)
        call CMTSVT(x_p,u0,const,xdot_p2,y_p2)
        A(:,k)=(xdot_p1-xdot_p2)/(2*DELXLIN(k))
        C(:,k)=(y_p1-y_p2)/(2*DELXLIN(k))
        
    end do
    
    ! control matrix
    do k = 1,NCTRLS
        u_p=u0
        u_p(k,1)=u_p(k,1)+DELCLIN(k)
        call CMTSVT(x0,u_p,const,xdot_p1,y_p1)
        u_p(k,1)=u_p(k,1)-2*DELCLIN(k)
        call CMTSVT(x0,u_p,const,xdot_p2,y_p2)
        B(:,k)=(xdot_p1-xdot_p2)/(2*DELCLIN(k));
        D(:,k)=(y_p1-y_p2)/(2*DELCLIN(k));
    end do

end subroutine