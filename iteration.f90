program iteration

  implicit none
  
  double precision, dimension (0:100) :: x
  integer :: n
  
  x(0) = 0.3d0
  
  call newtonRaphson
  
  contains
  
  subroutine pontoFixo
	
	implicit none

	double precision :: tol
	integer :: n
	
	tol = 1.0d-15
	
	n = 0
	
	do while (abs(g(x(n)) - x(n)) > tol .AND. n < size(x))
	  
	  write(*,*) n, x(n)
      x(n + 1) = g(x(n))
      n = n + 1
      
    end do
	
  end subroutine pontoFixo
  
  subroutine newtonRaphson
  
    implicit none
    
    double precision :: tol
	integer :: n
	
	tol = 1.0d-15
	
	n = 0
    
    do while (n < size(x))
    
       x(n + 1) = x(n) - f(x(n))/flinha(x(n))
       
       write(*,*) n, x(n)
       
       if (abs(x(n + 1) - x(n)) < tol) exit
       n = n + 1
       
    end do
  
  end subroutine
  
  function g(x)

    implicit none

    double precision :: g, x
    
    !g = cos(x)**(2.0d0/3.0d0)
    !g = (7*(x**2) + 14*x + 6)**(1.0/3.0)
    
    return
    
  end function g

  function f(x)

    implicit none

    double precision :: f,x
    
    f = (1.0/3.0)*(x**3.d0) - 2.d0*(x**2.d0) + 3.d0*x
    
    return
    
  end function f

  function flinha(x)

    implicit none

    double precision :: flinha, x 
    
    flinha = x**2.d0 - 4.d0*x + 3
    
    return
    
  end function flinha

end program iteration
