program derivada

  implicit none

  double precision :: h
  integer :: i, N
  double precision, dimension (:), allocatable :: x, flinhaav, flinhaat

  N = 40
  h = 1.d0/dble(N)

  allocate(x(0:N))
  allocate(flinhaav(0:N-1)) !para guardar os resultados da diferença finita AVançada
  allocate(flinhaat(1:N)) !para guardar os resultados da diferença finita ATrasada
  
  !Construção dos pontos da tabela no intervalo de x=0 a x=1
  do i=0,N
     x(i)=i*h
  end do

  open(unit=123,file='saida_av.dat', status='unknown')
  open(unit=124,file='saida_at.dat', status='unknown')

  !Cálculo das derivadas
  do i=0,N-1
     flinhaav(i) = (f(x(i+1))-f(x(i)))/h !diferença finita avançada
     write(123,*) x(i), flinhaav(i) 
  end do
  do i=1,N
     flinhaat(i) = (f(x(i))-f(x(i-1)))/h  !diferença finita atrasada
     write(124,*) x(i), flinhaat(i)
  end do

  !Esta linha é para análise da ordem. Faremos em sala juntos. 
  !write(*,*) h, x(N/2), flinhaav(N/2), flinhaat(N/2)

  close(unit=123)
  close(unit=124)

  deallocate(x)
  deallocate(flinhaav)
  deallocate(flinhaat)
  
contains

  function f(x)

    implicit none

    double precision :: f, x
    
    f = x*dsin(x)
    
    return
    
  end function f
  
end program derivada

