module func

   implicit none

   contains 
   ! Returns the inverse of a matrix calculated by finding the LU
   ! decomposition.  Depends on LAPACK.
   function inv(A) result(Ainv)
      real*8, dimension(:,:), intent(in) :: A
      real*8, dimension(size(A,1),size(A,2)) :: Ainv
   
      real*8, dimension(size(A,1)) :: work  ! work array for LAPACK
      integer, dimension(size(A,1)) :: ipiv   ! pivot indices
      integer :: n, info
   
      ! External procedures defined in LAPACK
      external DGETRF
      external DGETRI
   
      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      Ainv = A
      n = size(A,1)
   
      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call DGETRF(n, n, Ainv, n, ipiv, info)
   
      if (info /= 0) then
         stop 'Matrix is numerically singular!'
      end if
   
      ! DGETRI computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
      call DGETRI(n, Ainv, n, ipiv, work, n, info)
   
      if (info /= 0) then
         stop 'Matrix inversion failed!'
      end if
   end function inv

   !!!!!! ------------------- this is hard coded ------------------- !!!!!

   function squeeze(A) result(B)
      real*8 :: A(3,3,2,2,10)
      real*8 :: B(3,3,10)
      !integer :: siz(size(A))
      integer :: i
      !logical :: is_scalar

      do i = 1,10
         B(:,:,i) = A(:,:,2,1,i)
      end do

      !siz = shape(A)
      !do i = 1, size(siz)
      !   if (siz(i) == 1) then
      !      siz(i) = 0
      !   end if
      !end do

      !B = reshape(A, siz)

      !A(:,:,2) = 0

      !siz_1D = size(A,Dim=1)
      !siz_2D = size(A,Dim=2)
      !siz_3D = size(A,Dim=3)
  
      !siz = (/siz_1D,siz_2D,siz_3D/)
  
      !do i = 1, size(A)
      !    if (sum(A(:,:,i)) == 0) then
      !        siz(i) = 0
      !    end if 
      !    if (siz(i) /= 0) then
      !        new_siz = siz(i)
      !    end if 
      !end do

   end 

   function permute(A,order) result(per_A) 
      real*8 :: A(3,3,10)
      integer :: order(3)
      real*8 :: per_A(10,3,3)
      !real*8, allocatable :: B(:,:,:)
      
      !integer :: siz_1D, siz_2D, siz_3D
      integer :: i,j,k
  
  
      do i = 1,3
          do j = 1,3
              do k = 1,10
                  per_A(k,i,j) = A(i,j,k)
              end do
          end do
      end do
  
  
      !siz_1D = size(A,Dim=1)
      !siz_2D = size(A,Dim=2)
      !siz_3D = size(A,Dim=3)
  
      !allocate(B())
   end function

   function interp1(x,v,xq) result(vout)

      ! linear interpolation for single query point and 3D v
      real*8 :: x(10), v(10,3,3), xq
      real*8 :: vout(1,3,3)
      integer :: i,j,k
  
      do k = 1,size(v,dim=3)
          do j = 1,size(v,dim=2)
              do i = 1,size(x)
                  if (x(i) == xq) then
                      vout(1,j,k) = v(i,j,k)
                  end if
              end do
          end do
      end do
  
  
  end function

end module func