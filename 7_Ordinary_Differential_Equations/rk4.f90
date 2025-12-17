program rk4

  implicit none

  double precision, dimension (:), allocatable :: y
  double precision :: k1, k2, k3, k4, h
  integer :: n, Nfim, k

  open(unit=124,file='ordem.dat', status='unknown')

  do k = 0, 10

     Nfim = 20*(2**k)
     h = 2.d0/dble(Nfim)

     allocate(y(0:Nfim))

     y(0) = 1.d0
     do n = 0,Nfim-1
        k1 = h*f(n*h,y(n))
        k2 = h*f(n*h+0.5d0*h,y(n)+0.5d0*k1)
        k3 = h*f(n*h+0.5d0*h,y(n)+0.5d0*k2)
        k4 = h*f(n*h+h,y(n)+k3)
        y(n+1) = y(n) + (k1+2.d0*k2+2.d0*k3+k4)/6.d0
     end do

!!$  open(unit=123,file='saida.dat',status='unknown')
!!$  do n=0,Nfim
!!$     write(123,*) n*h,y(n)
!!$  end do
!!$  close(unit=123)  

     write(124,*) h, dabs(y(Nfim)-dexp(-2.d0/3.d0))


     deallocate(y)

  end do

  close(unit=124)
contains

  function f(t,y)

    implicit none

    double precision :: f, y, t

    f = (1.d0-4.d0*t/3.d0)*y

    return

  end function f


end program rk4

