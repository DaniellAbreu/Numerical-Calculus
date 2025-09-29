program bissec

  implicit none

  double precision :: xmed, a, b, tol, erro
  integer :: n
  

  a = 0.d0
  b = 6.d0
  tol = 1.d-15
  n = 0
  erro = 1233.d0
  
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
    
    f = x**3.d0 - x - 1.d0
    
    return
    
  end function f
  
end program bissec

