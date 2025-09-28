program matrix
 
  implicit none
 
  integer :: i,j,k,n
  double precision, dimension (:,:), allocatable :: m, lambda
 
  n = 5
  allocate(lambda(n,n))
  allocate(m(n,n+1))
 
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
 
  write(*,*) 'Depois'
  do i=1,n
     write(*,*) (m(i,j),j=1,n+1)
  end do
 
  
  deallocate(lambda)
  deallocate(m)
 
end program matrix
