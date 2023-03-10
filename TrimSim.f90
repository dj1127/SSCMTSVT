subroutine TrimSim(aircraft,x0,u0,targ_des,XSCALE,YSCALE,TRIMVARS,&
    TRIMTARG,NSTATES,NCTRLS,x0trim,u0trim,itrim)

    use defineAircraftProperties

    implicit none

    type(const_struct) :: const

    character(16), intent(in) :: aircraft
    real*8,intent(in) :: x0(12,1),u0(20,1),targ_des(14,1),
    real*8,intent(out) :: x0trim,u0trim,itrim
    real*8 :: XSCALE(12,1),YSCALE(10,1),TRIMVARS(14),&
    TRIMTARG(14),NSTATES,NCTRLS


    real*8 :: error, targ_err, trim_tol
    real*8 :: XYSCALE(12,1)
    integer :: it, itmax

    ! ---------CMTSVT---------
    real*8 :: y0, xdot0
    real*8, allocatable:: temp(:)


    ! Unpack aircraft constants
    XSCALE=const%XSCALE;
    YSCALE=const%YSCALE; 
    TRIMVARS=const%TRIMVARS;
    TRIMTARG=const%TRIMTARG;
    NSTATES=const%NSTATES;
    NCTRLS=const%NCTRLS;

    ! initial state and control input guess
    x0trim = x0
    u0trim = u0
    ! initialize number of iterations and error
    it = 0
    error = 100
    ! tolerance on trim error
    trim_tol = 5e-4
    ! maximum number of iterations
    itmax = 100
    ! trim aircraft
    write(*, '(A)') 'ITERATION      TRIM ERROR'
    do while ((it < itmax) .and. (error > trim_tol))
        it = it + 1
        call CMTSVT(x0trim, u0trim, const, xdot0, y)
        temp = [xdot0, y0]
        targvec = temp(TRIMTARG)
        targ_err = targvec - targ_des
        XYSCALE = [XSCALE, YSCALE]
        error = maxval(abs(targ_err) / XYSCALE(TRIMTARG))
        write(*, '(I2, 2X, F5.4)') it, error
        if (err > trim_tol) then
            call LinSim(A, B, C, D, aircraft, x0trim, u0trim, const)
            Jac = [A, B; C, D]
            Jac = Jac(TRIMTARG, TRIMVARS)
            trimvec = [x0trim, u0trim]
            trimvec(TRIMVARS) = trimvec(TRIMVARS) - 0.5 * matmul(pinv(Jac), targ_err)
            x0trim = trimvec(1:NSTATES)
            u0trim = trimvec(NSTATES+1:NSTATES+NCTRLS)
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