program raizes

  implicit none

  double precision, dimension (0:100) :: x
  integer :: n
  
  !x(0) = 1.d0
  x(0) = 3.d0
  do n=0,99
     !x(n+1) = (x(n)+1.0d0)**(1.0d0/3.0d0)
     !x(n+1) = x(n)**3-1.0
     x(n+1) = g(x(n))
  end do

  do n=0,100
     write(*,*) n, x(n)
  end do

contains

  function g(x)

    implicit none

    double precision :: g, x
    
    !g = (x+1.d0)**(1.d0/3.d0)

    g = x - f(x)/flinha(x)
    
    return
    
  end function g


  function f(x)

    implicit none

    double precision :: f,x
    
    f = x**3.d0 - x - 1.d0
    
    return
    
  end function f

  function flinha(x)

    implicit none

    double precision :: flinha, x 
    
    flinha = 3.d0*x**2.d0 - 1.d0
    
    return
    
  end function flinha
  
end program raizes

