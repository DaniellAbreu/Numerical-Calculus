program iterativo

  implicit none

  double precision :: x
  integer :: i, N

  N = 100

  open(unit=123,file='saida.dat',status='unknown')
  x =34.d0
  do i=1,N
     x = g(x)
     write(123,*) i, x
  end do
  close(unit=123)
  
contains

  function g(x)

    implicit none

    double precision :: g, x

    !g = dcos(x)
    g = x*x
    
    return
    
  end function g
  

end program iterativo

