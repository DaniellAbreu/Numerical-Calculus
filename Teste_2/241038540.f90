! Compilador: gfortran
! Compilação: gfortran 241038540.f90
! Execução: ./a.out

program teste_2

  implicit none
  
  double precision, dimension (0:100) :: x
  integer :: n
  
  x(0) = 0.7d0 ! Chute Inicial.
  
  call newtonRaphson
  
  contains
  
  subroutine newtonRaphson
  
    implicit none
    
    double precision :: tol
	integer :: n
	
	tol = 1.0d-15 ! Tolerância do Erro.
	
	n = 0
    
    do while (n < size(x))
    
       x(n + 1) = x(n) - f(x(n))/flinha(x(n))
       
       write(*,*) n, x(n) ! Imprimindo Número de Iterações até Valor de Erro Menor que a Tolerância.
       
       ! Última Iteração será a Raiz da Função dentro do Intervalo e da Tolerância.
       
       if (abs(x(n + 1) - x(n)) < tol) exit
       n = n + 1
       
    end do
  
  end subroutine

  function f(x)

    implicit none

    double precision :: f,x
    
    f = cos(x) - acos(-1.d0)/x
    
    return
    
  end function f

  function flinha(x)

    implicit none

    double precision :: flinha, x 
    
    flinha = acos(-1.d0)*(x**(-2.d0)) - sin(x)
    
    return
    
  end function flinha

end program teste_2
