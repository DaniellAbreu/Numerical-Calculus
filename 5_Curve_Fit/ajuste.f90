program ajuste

  implicit none


  double precision, dimension (:,:), allocatable:: vand, vandtransp, gram
  double precision, dimension (:), allocatable :: t, y, atb, coefs
  double precision :: soma
  integer :: i, j, k, grau, N

  N = 5
  grau = 2 
  
  allocate(vand(N+1,grau+1))
  allocate(vandtransp(grau+1,N+1))
  allocate(gram(grau+1,grau+1))
  allocate(y(N+1))
  allocate(t(N+1))
  allocate(atb(grau+1))
  allocate(coefs(grau+1))

 
  !leitura dos dados
  !open(unit=123,file='dados.dat', status='old')
  !do i = 1,N+1
  !   read(123,*) t(i), y(i)
  !   !write(*,*) i, t(i), y(i)
  !end do
  !close(unit=123)

  t(1) = 1.0d0
  t(2) = 2.0d0
  t(3) = 3.0d0
  t(4) = 4.0d0
  t(5) = 5.0d0
  t(6) = 6.0d0
  
  y(1) = 0.987d0
  y(2) = 1.46d0
  y(3) = 3.18d0
  y(4) = 3.18d0
  y(5) = 6.7d0
  y(6) = 4.8d0

  
  !Gerar Vandermonde
  do i=1,N+1
     do j=1,grau+1
        vand(i,j) = t(i)**(j-1)
     end do
  end do

  !cáculo da transposta
  do i=1,N+1
     do j=1,grau+1
        vandtransp(j,i) = vand(i,j)
     end do
  end do

  !multiplicação de matrizes para **regularização do sistema**
  do i=1,grau+1
     do j = 1, grau+1
        soma = 0.d0
        do k=1,N+1
           soma = soma + vandtransp(i,k)*vand(k,j)
        end do
        gram(i,j) = soma
     end do
  end do

  do i=1,grau+1
     write(*,*) (gram(i,j),j=1,grau+1)
  end do

  !multiplicação da transposta pelo vetor solução para **regularização do sistema**
  do i=1,grau+1
     soma = 0.d0
     do k=1,N+1
        soma = soma + vandtransp(i,k)*y(k)
        !write(*,*) i,k,soma, y(k)
     end do
     atb(i) = soma
  end do
  
  !cálculo dos mínimos quadrados -- AJUSTE
  call sistema(gram,atb,coefs,grau)

  !cálculo dos polinômio interpolador
  !call sistema(vand,y,coefs,grau)

  !impressão em tela do resultado -- coeficientes do polinômio
  write(*,*) coefs

  deallocate(vand)
  deallocate(vandtransp)
  deallocate(gram)
  deallocate(y)
  deallocate(t)
  deallocate(atb)
  deallocate(coefs)
  
contains
  
  subroutine sistema(a,b,x,grau)

    implicit none

    double precision, dimension(grau+1,grau+1) :: a
    double precision, dimension(grau+1) :: x,b
    double precision :: soma
    integer :: grau, i,j,k
    
    do k = 1, 40000
       do i=1, grau+1
          soma = 0.d0
          do j=1,grau+1
             if(j.ne.i) then
                soma = soma + a(i,j)*x(j)
             end if
          end do
          x(i) = (b(i) - soma)/a(i,i)
       end do
    end do

    return
    
  end subroutine sistema

  
end program ajuste

