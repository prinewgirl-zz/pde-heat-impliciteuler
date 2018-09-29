program implicit_euler
      implicit none    
      real(kind=8) :: alpha,t0,tf            
      integer :: i,j, times
      integer n, nrhs,ldb,info
      double precision, allocatable :: dl(:),d(:),du(:), b(:),t(:)
      
            
      write(*,*) ' type  n, lambda, times units (times), t0 e tf'
      read (*,*) n,alpha, times, t0, tf
      
     allocate (dl(n-1))
     allocate (d(n))
     allocate (du(n-1))  
     allocate (b(n))
     allocate (t(n+2))
     
     nrhs = 1
     ldb = n
     

   !Superdiagonal,
   !Diagonal,
   !Subdiagonal.
    du(1:n-1) = -alpha
    d(1:n) = 1.0+2*alpha
    dl(1:n-1) = -alpha
    !initial conditions    
     t(0) = t0
     t(1:n-1) = 0.0   
     t(n+1) = tf
    
     do j = 1,times     
       b(1) = t(1)+alpha*t(0)
       do i = 2, n-1
         b(i) = t(i)
       end do     
       b(n) = t(n)+alpha*t(n+1)
   
       call dgtsv(n,nrhs,dl,d,du,b,ldb,info)
       write(*,*) b
       t(1:n) = b
       
      end do
     
    
endprogram implicit_euler
