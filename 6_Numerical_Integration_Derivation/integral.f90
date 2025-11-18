program integral

  implicit none

  double precision :: a, b, h, soma
  double precision :: t1,t2,t3,a1,a2,a3
  double precision, dimension (:), allocatable :: t
  integer :: i,N
  
  N = 20

  allocate(t(0:N))

  !vamos calcular a integral de x^5 entre 0 e 1.
  a = 0.d0 !limites de integração
  b = 1.d0 !limites de integração
  h = dabs(b-a)/dble(N) !passo de integração

  !Construção dos pontos onde a quadratura vai ser avaliada
  do i=0,N
     t(i) = a + i*h
  end do

  !Retângulos pela esquerda
  soma = 0.d0
  do i=0,N-1
     soma = soma + f(t(i))*h
  end do
  write(*,*) 'A integral é (retangulo)', soma
  
  !Regra do Trapézio
  soma = 0.d0
  do i=1,N-1
     soma = soma + f(t(i))
  end do
  soma = h*(soma + (f(t(0))+f(t(N)))/2.d0)
  write(*,*) 'A integral é (trapezio)', soma

  !Nós da quadratura gaussiana para P3(t)) - Retirados de wikipedia 
  !https://en.wikipedia.org/wiki/Gaussian_quadrature
  !Atenção: aqui já foi feita a translação do domínio de integração
  t1 = ((b-a)/2.d0)*(-dsqrt(3.d0/5.d0))+((b+a)/2.d0)
  t2 = ((b-a)/2.d0)*0.d0+((b+a)/2.d0)
  t3 = ((b-a)/2.d0)*(dsqrt(3.d0/5.d0))+((b+a)/2.d0)

  !Pesos da quadratura gaussiana para P3(t) - Retirados de wikipedia 
  !https://en.wikipedia.org/wiki/Gaussian_quadrature
  !Vocês devem saber como calcular!
  a1 = 5.d0/9.d0
  a2 = 8.d0/9.d0
  a3 = 5.d0/9.d0

  !Cálculo da quadratura
  soma = (a1*f(t1)+a2*f(t2)+a3*f(t3))*((b-a)/2.d0)
  write(*,*) 'A integral é (gaussiana)', soma

  write(*,*) 'A integral é (exata)', 1.d0/6.d0
  
  deallocate(t)
  
contains

  function f(x)

    implicit none

    double precision :: x, f
   
    f = x**5.0
    
    return
    
  end function f
  
end program integral

