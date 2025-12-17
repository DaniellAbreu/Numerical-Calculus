program euler

  implicit none

  double precision, dimension (:), allocatable :: y
  integer :: n, Nfim, k
  double precision :: h
  
  open(unit=124,file='ordem.dat',status='unknown')

  do k = 0, 10
     Nfim = 20*(2**k)
     h = 2.d0/dble(Nfim)

     allocate(y(0:Nfim))

     y(0) = 1.d0
     do n = 0, Nfim-1
        y(n+1) = y(n) + h*f(n*h,y(n)) !t(n) = n*h
     end do

     write(124,*) h, dabs(y(Nfim)-dexp(-2.d0/3.d0))

!!$     open(unit=123,file='saida.dat',status='unknown')
!!$     do n=0,Nfim
!!$        write(123,*) n*h,y(n)
!!$     end do
!!$     close(unit=123)

     deallocate(y)
  end do

  close(unit=124)
  
   contains

     function f(t,y)

       implicit none


       double precision :: f, t, y

       f = (1.d0 - 4.d0*t/3.d0)*y

       return

     end function f

   end program euler

