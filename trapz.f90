! Trapezodial Integration

    subroutine trapz(rng,fn,d,res)
        IMPLICIT none
        integer :: d
        real*8 :: rng(d),fn(d)
        real*8 :: rng_t(d,1),fn_t(d,1)
        real*8 :: diff(d-1,1), fn_calc1(d-1,1),fn_calc2(d-1,1)
        real*8 :: diff_t(1,d-1)
        real*8 :: res
        integer :: nshifts,i
        !integer, dimension (d:1) :: order 
        
        ! initialize res
        res = 0.
        
        nshifts = 0
        
        rng_t = reshape(rng,(/ d, 1 /))
        fn_t = reshape(fn,(/ d, 1 /))

        do i = 1,d-1
            diff = rng(i+1)-rng(i)
        end do

        diff_t = transpose(diff)
        
        do i =1,d-1
            fn_calc1(i,1) = fn_t(i,1)
        end do
        
        do i =2,d
            fn_calc2(i-1,1) = fn_t(i,1)
            
        end do
        
        do i = 1,d-1
            res = res + diff_t(1,i)*(fn_calc1(i,1)+fn_calc2(i,1))/2

        end do

    end
    
