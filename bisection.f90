program bissec

  implicit none

  double precision :: xmed, a, b, tol, erro
  integer :: n
  

  a = 1.d0
  b = 6.d0
  tol = 1.d-8
  n = 0
  erro = 1234.d0
  
  do while (erro > tol)
     n = n + 1
     xmed = (a+b)/2.d0 !ponto médio
     erro = dabs(f(xmed))
     if(erro > tol) then
        if(f(a)*f(xmed) < 0.d0) then
           b = xmed
        else
           a = xmed
        end if
     end if
     write(*,*) n, xmed, a, b
  end do
  
contains

  function f(x)

    implicit none

    double precision :: f, x
    
    f = x**4.d0 - 23.d0*(x**3.d0) + 179.d0*(x**2.d0) - 517.d0*x + 360.d0
    
    return
    
  end function f
  
end program bissec

