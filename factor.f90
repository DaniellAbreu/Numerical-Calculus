program factor

  implicit none

  integer :: i,j,k,n
  double precision, dimension (:,:), allocatable :: m, lambda, mold, u, inversa
  double precision, dimension (:), allocatable :: x, b, c
  double precision :: soma
  
  n = 5
  
  allocate(lambda(n,n))
  allocate(u(n,n))
  allocate(inversa(n,n))
  allocate(m(n,n+1))
  allocate(mold(n,n+1))
  allocate(x(n))
  allocate(b(n))
  allocate(c(n))
  

  m(1,1) = 2.d0
  m(1,2) = 1.d0
  m(1,3) = -1.d0
  m(1,4) = 1.d0
  m(1,5) = 3.d0
  m(1,6) = 7.d0

  m(2,1) = 1.d0
  m(2,2) = 0.d0
  m(2,3) = 2.d0
  m(2,4) = -1.d0
  m(2,5) = 1.d0
  m(2,6) = 2.d0

  m(3,1) = 0.d0
  m(3,2) = -2.d0
  m(3,3) = 1.d0
  m(3,4) = 1.d0
  m(3,5) = -1.d0
  m(3,6) =  2.d0

  m(4,1) = 3.d0
  m(4,2) = 1.d0
  m(4,3) = -4.d0
  m(4,4) = 0.d0
  m(4,5) = 5.d0
  m(4,6) = 6.d0

  m(5,1) = 1.d0 
  m(5,2) = -1.d0 
  m(5,3) = -1.d0
  m(5,4) = -1.d0
  m(5,5) = 1.d0
  m(5,6) = 3.d0

  mold = m

  lambda = 0.d0
  do i=1,n
     lambda(i,i) = 1.d0
  end do
  
  write(*,*) 'Antes'
  do i=1,n
     write(*,*) (m(i,j),j=1,n+1) !'do' implícito
  end do
  
  do j = 1, n !percorremos as n primeiras colunas da matriz
     if(m(j,j).eq.0.d0) then
        write(*,*) 'Problema! Pivô nulo linha', j
        stop
     else
        do i=j+1,n !percorremos as linhas abaixo do pivô
           lambda(i,j) = m(i,j)/m(j,j)
           do k=1,n+1 !lembrar da matriz extendida/combinação linear
              m(i,k) = m(i,k) - lambda(i,j)*m(j,k)
           end do
        end do
     end if
  end do

  !construção da matriz u
  do i=1,n
     do j=1,n
        u(i,j) = m(i,j)
     end do
  end do

  write(*,*) 'Depois'
  do i=1,n
     write(*,*) (m(i,j),j=1,n+1)
  end do

  !algoritmo de substituição reversa
  x(n) = m(n,n+1)/m(n,n)
  do i = n-1, 1, -1
     soma = 0.d0
     do k = i+1,n !soma dos termos à direita da diagonal
        soma = soma + m(i,k)*x(k)
     end do
     x(i) = (m(i,n+1)-soma)/m(i,i)
  end do

  write(*,*) 'Solução'
  write(*,*) x

  !rotina de verificação
  do i = 1, n
     soma = 0.d0
     do j=1,n
        soma = soma + mold(i,j)*x(j)
     end do
     write(*,*) 'linha', i, dabs(mold(i,n+1)-soma)
  end do

  write(*,*) 'Matriz U'
  do i=1,n
     write(*,*) (u(i,j),j=1,n)
  end do
  write(*,*) 'Matriz L'
  do i=1,n
     write(*,*) (lambda(i,j),j=1,n)
  end do

  write(*,*) 'verificação LU'
  do i=1,n
     do j=1,n
        soma = 0.d0
        do k=1,n
           soma = soma + lambda(i,k)*u(k,j)
        end do
        write(*,*) i,j, dabs(mold(i,j)-soma)
     end do
  end do

  !cálculo da inversa de A

  do j=1,n !para cada uma das colunas da inversa
     b = 0.d0
     b(j) = 1.d0
     !algoritmo de substituição direta
     c(1) = b(1)
     do i=2,n
        soma = 0.d0
        do k=1,i-1
           soma = soma + lambda(i,k)*c(k)
        end do
        c(i) = b(i) - soma
     end do
     
     !algoritmo de substituição reversa
     x(n) = c(n)/u(n,n)
     do i = n-1, 1, -1
        soma = 0.d0
        do k = i+1,n !soma dos termos à direita da diagonal
           soma = soma + u(i,k)*x(k)
        end do
        x(i) = (c(i)-soma)/u(i,i)
     end do
     do k = 1, n
        inversa(k,j) = x(k)
     end do
  end do

  write(*,*) 'Inversa'
  do i=1,n
     write(*,*) (inversa(i,j),j=1,n)
  end do

  write(*,*) 'verificação inversa'
  do i=1,n
     do j=1,n
        soma = 0.d0
        do k=1,n
           soma = soma + inversa(i,k)*mold(k,j)
        end do
        write(*,*) i,j, soma
     end do
  end do
  
  deallocate(lambda)
  deallocate(m)
  deallocate(mold)
  deallocate(x)
  deallocate(u)
  deallocate(inversa)
  deallocate(b)
  deallocate(c)
  
end program factor
