subroutine TrimSim(aircraft,x0,u0,targ_des,const,x0trim,u0trim,itrim)

    use AllType

    implicit none

    type(const_struct) :: const

    character(16), intent(in) :: aircraft
    real*8,intent(in) :: x0(12,1),u0(20,1),targ_des(14,1)
    real*8,intent(out) :: x0trim(12,1),u0trim(20,1),itrim
    real*8,allocatable :: XSCALE(:,:),YSCALE(:,:),TRIMVARS(:),TRIMTARG(:),Jac_temp(:,:),Jac(:,:),Jac_mul(:,:),Jac_pinv(:,:)
    real*8 :: NSTATES,NCTRLS
    real*8 :: A(12,12),B(12,20),C(10,12),D(10,20)



    real*8 :: err, trim_tol
    real*8 :: XYSCALE(12,1)
    integer :: it, itmax

    ! ---------CMTSVT---------
    real*8 :: y0(10,1), xdot0(12,1)
    real*8,allocatable :: temp(:,:), targvec(:,:), targ_err(:,:)
    real*8,allocatable :: trimvec(:,:)

    ! External Subroutines
    EXTERNAL DGETRF


    

    ! Unpack aircraft constants
    XSCALE=const%XSCALE
    YSCALE=const%YSCALE
    TRIMVARS=const%TRIMVARS;
    TRIMTARG=const%TRIMTARG;
    NSTATES=const%NSTATES;
    NCTRLS=const%NCTRLS;

    ! allocate array size
    allocate(XSCALE(size(const%XSCALE,DIM=1),size(const%XSCALE,DIM=2)))
    allocate(YSCALE(size(const%YSCALE,DIM=1),size(const%YSCALE,DIM=2)))
    allocate(TRIMVARS(size(const%TRIMVARS)))
    allocate(TRIMTARG(size(const%TRIMTARG)))
    allocate(temp(size(xdot0, DIM = 1)+size(y0, DIM = 1),1))
    allocate(targvec(size(temp,DIM=1),1))
    allocate(targ_err(size(targvec,DIM=1),1))
    allocate(trimvec(size(x0, DIM = 1)+size(u0, DIM = 1),1))
    allocate(Jac_temp((size(A,Dim=1)+size(C,Dim=1)), (size(A,Dim=2)+size(B,Dim=2))))
    allocate(Jac((size(A,Dim=1)+size(C,Dim=1)), (size(A,Dim=2)+size(B,Dim=2))))
    allocate(Jac_mul((size(A,Dim=1)+size(C,Dim=1)), (size(A,Dim=2)+size(B,Dim=2))))
    allocate(Jac_pinv((size(A,Dim=1)+size(C,Dim=1)), 1))

    ! initial state and control input guess
    x0trim = x0
    u0trim = u0
    ! initialize number of iterations and error
    it = 0
    err = 100
    ! tolerance on trim error
    trim_tol = 5e-4
    ! maximum number of iterations
    itmax = 100
    ! trim aircraft
    write(*, '(A)') 'ITERATION      TRIM ERROR'
    do while ((it < itmax) .and. (err > trim_tol))
        it = it + 1
        call CMTSVT(x0trim, u0trim, const, xdot0, y0)
        temp = reshape((/xdot0, y0/),(/22,1/)) ! need to be edited
        targvec = temp(TRIMTARG,:)
        targ_err = targvec - targ_des
        XYSCALE = reshape((/XSCALE, YSCALE/),(/12,1/))
        err = maxval(abs(targ_err(:,1)) / XYSCALE(TRIMTARG,1))
        write(*, '(I2, 2X, F5.4)') it, err
        if (err > trim_tol) then
            call LinSim(x0trim, u0trim, const, A, B, C, D)
            Jac_temp = reshape((/A, B, C, D/),(/ (size(A,Dim=1)+size(C,Dim=1)), (size(A,Dim=2)+size(B,Dim=2)) /),&
             order = (/ 2, 1 /))
            Jac = Jac_temp(TRIMTARG, TRIMVARS)
            trimvec = reshape((/x0trim, u0trim/),(/size(x0,DIM=1)+size(u0,DIM=1),1/))
            ! pseudo inverse
            Jac_mul = matmul(transpose(Jac),Jac)
            call DGETRF(size(A,Dim=1),size(A,Dim=2),Jac_mul,size(A,Dim=1))
            Jac_pinv = matmul(Jac_mul,transpose(Jac))

            trimvec(TRIMVARS,:) = trimvec(TRIMVARS,:) - 0.5 * matmul(Jac_pinv, targ_err)
            x0trim = trimvec(1:NSTATES,:)
            u0trim = trimvec(NSTATES+1:NSTATES+NCTRLS,:)
        end if
    end do
    ! store info on whether trim was achieved or not 
    if (err > trim_tol) then
        write(*, '(A)') 'Warning: trim not achieved'
        itrim = 0
    else
        write(*, '(A)') 'Successful trim'
        itrim = 1
    end if

end