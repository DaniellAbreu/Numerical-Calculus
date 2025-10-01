program sist_iter

  implicit none

  integer :: i, j, k, N
  double precision :: soma, tol
  double precision, dimension (:,:), allocatable :: a
  double precision, dimension (:), allocatable :: b, x, xold

  N = 5
  tol = 1.d-14
  
  allocate(a(N,N))
  allocate(b(N))
  allocate(x(N))
  allocate(xold(N))

!!$  a(1,1) = 1.d0
!!$  a(1,2) = 2.d0
!!$  a(1,3) = 1.d0
!!$
!!$  a(2,1) = 2.d0
!!$  a(2,2) = 6.d0
!!$  a(2,3) = 1.d0
!!$
!!$  a(3,1) = 1.d0
!!$  a(3,2) = 1.d0
!!$  a(3,3) = 4.d0
!!$
!!$  b(1) = 2.d0
!!$  b(2) = 7.d0
!!$  b(3) = 3.d0

  a(1,1) = 8.d0
  a(1,2) = 1.d0
  a(1,3) = -1.d0
  a(1,4) = 1.d0
  a(1,5) = 3.d0
  b(1) = 7.d0

  a(4,1) = 1.d0
  a(4,2) = 0.d0
  a(4,3) = 2.d0
  a(4,4) = -8.d0
  a(4,5) = 1.d0
  b(4) = 2.d0

  a(3,1) = 0.d0
  a(3,2) = -2.d0
  a(3,3) = 7.d0
  a(3,4) = 1.d0
  a(3,5) = -1.d0
  b(3) =  2.d0

  a(2,1) = 3.d0
  a(2,2) = 18.d0
  a(2,3) = -4.d0
  a(2,4) = 0.d0
  a(2,5) = 5.d0
  b(2) = 6.d0

  a(5,1) = 1.d0 
  a(5,2) = -1.d0 
  a(5,3) = -1.d0
  a(5,4) = -1.d0
  a(5,5) = 9.d0
  b(5) = 3.d0

  x(1) = 0.d0 !Chute inicial
  x(2) = 0.d0 !Chute inicial
  x(3) = 0.d0 !Chute inicial
  x(4) = 0.d0 !Chute inicial
  x(5) = 0.d0 !Chute inicial 
  
  !Notação: x = (x(1), x(2), x(3)) = (x, y, z)

  k = 0
  do while (residuo(b,a,x,N) .ge. tol)
     k = k + 1
     !xold = x !Cópia de x para Gauss-Jacobi
     do i=1,N !percorrer todas as linhas do sistema
        soma = 0.d0
        do j=1,N !percorrer todas as colunas do sistema
           if(j.ne.i) then !apenas termos não-diagonais
              soma = soma + a(i,j)*x(j) !xold para Gauss-Jacobi, x para Gauss-Seidel
           end if
        end do
        x(i) = (b(i)-soma)/a(i,i)
     end do
     write(*,*) k, residuo(b,a,x,N)
  end do
  write(*,*) 'A solução é', x


  deallocate(a)
  deallocate(b)
  deallocate(x)
  deallocate(xold)

contains

  function residuo(b,a,x,N)

    implicit none

    integer :: i, N, j
    double precision, dimension (N,N) :: a
    double precision, dimension (N) :: b,x,c
    double precision :: residuo, soma
    
    do i = 1, N
       soma = 0.d0
       do j=1,N
          soma = soma + a(i,j)*x(j)
       end do
       c(i) = b(i) - soma
    end do

    residuo = 0.d0
    do i = 1, N
       residuo = residuo + c(i)*c(i)
    end do
    residuo = dsqrt(residuo)
    
    return
    
  end function residuo
  
end program sist_iter

