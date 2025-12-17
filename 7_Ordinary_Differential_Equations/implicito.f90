program implicito

  implicit none

  double precision, dimension (:), allocatable :: y
  double precision :: h
  integer :: n, Nfim, k
  
  Nfim = 80
  h = 2.d0/dble(Nfim)

  allocate(y(0:Nfim))

  open(unit=123, file='saida.dat',status='unknown')
  y(0) = 1.d0
  do n = 0, Nfim-1 ! laço de Euler
     do k = 1, 10 !laço do Newton Raphson
        y(n+1) = y(n+1) - g(y(n+1),y(n),(n+1)*h, h)/glinha(y(n+1),y(n),(n+1)*h, h)
     end do
     write(123,*) (n+1)*h,y(n+1)
  end do
  close(unit=123)

  deallocate(y)
  
contains

  function g(y,yold,t,h)

    implicit none

    double precision :: y, yold, t, h, g

    g = y - yold - h*f(t,y)
    
    return
    
  end function g

  function glinha(y,yold,t,h)

    implicit none

    double precision :: y, yold, t, h, glinha

    glinha = 1.d0 - h*3.d0*(y**2.d0)*(1.d0-4.d0*t/3.d0)
    
    return
    
  end function glinha

  function f(t,y)

    implicit none

    double precision :: f, t, y 
    
    f = (1.d0-4.d0*t/3.d0)*(y**3.d0)

    return
    
  end function f


end program implicito

