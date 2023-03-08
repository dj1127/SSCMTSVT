! ---------------------------- TRIM PARAMETER -----------------------------

    module trim_parameter
    

    
    real*8, allocatable :: TRIMVARS(:), TRIMTARG(:)
    
    real*8 :: swashplate0(3), vel0(3), ang0(3),lambda0(3),lambdatot0(3)
    
    real*8 :: u0(20,1), x0(12,1)
    
    real*8 :: x0trim(12), y0trim(10), u0trim(20)
    
    real*8  :: x0trim_val(12), y0trim_val(10), u0trim_val(20)
    
    real*8 :: CTdes(2), Tdes(2), targ_des(14,1)
        
    end