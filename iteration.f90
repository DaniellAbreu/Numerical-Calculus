program iteration
  
  ! Métodos de Convergência de Pontos Fixos - Ponto Fixo, Newton-Raphson e Bisseção.
    
  implicit none
  
  double precision, dimension (:), allocatable :: x
  integer :: N
  
  N = 100
  
  allocate(x(0:N))
  
  ! Chute Inicial
  x(0) = 1.d0
  
  call newtonRaphson
  call pontoFixo
  call bissec
  
  deallocate(x)
  
  contains
  
  subroutine pontoFixo
	
	implicit none
	
	double precision :: tol
	integer :: n
	
	tol = 1.d-8
	
	n = 0
	
    do while (n < size(x) - 1)
    
      x(n + 1) = g(x(n))
      write(*,*) "Iteração Ponto Fixo:", n, x(n + 1)
      
      if (abs(x(n + 1) - x(n)) < tol) then
        write(*,*) "Convergiu para:", x(n+1)
        exit
      end if
      
      n = n + 1
    end do
	
  end subroutine pontoFixo
  
  subroutine newtonRaphson
  
    implicit none
    
    double precision :: tol
	integer :: n
	
	tol = 1.d-8
	
	n = 0
    
    do while (n < size(x) - 1)
    
      x(n + 1) = x(n) - f(x(n))/flinha(x(n))
       
      write(*,*) "Iteração Newton Raphson:", n, x(n + 1)
       
      if (abs(x(n + 1) - x(n)) < tol) then
        write(*,*) "Convergiu para:", x(n + 1)
        exit
      end if
       
      n = n + 1
    end do
  
  end subroutine
  
  subroutine bissec
  
    implicit none
    
    double precision :: tol, erro, xmed, a, b, fmid, fa
    integer :: n
    
    ! Intervalo: a > b
    a = 0.d0
    b = 6.d0
    
    n = 0
    
    tol = 1.d-8
    erro = 1000.d0
    
    fa = f(a)
    
    do while (erro > tol)

      xmed = (a + b)/2.d0
      fmid = f(xmed)
      
      write(*,*) "Iteração Bisseção", n, xmed
      
      n = n + 1
	  
	  if (fmid*fa < 0.d0) then
	    b = xmed
	  else
		a = xmed
	    fa = fmid
      end if
		
	  erro = dabs(fmid)
	  
    end do
    
    write(*,*) "Convergiu para:", xmed
  
  end subroutine
  
  function g(x)

    implicit none

    double precision :: g, x
    
    g = (x + 1)**(1.d0/3.d0)

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
    
    flinha = 3.d0*(x**2.d0) - 1.d0
    
    return
    
  end function flinha

end program iteration
