    subroutine ones(X,Y,matrix)
    
    integer :: X,Y
    real,dimension(X,Y) :: matrix
    
    do i=1,X
      do j = 1,Y
         matrix(i,j) = 1
      end do
    end do

    
    END subroutine