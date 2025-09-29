program vetores

  implicit none

  integer :: i, j, N
  real, dimension(:), allocatable :: x
  real, dimension(:,:), allocatable :: a, atransp

  N = 3

  allocate(x(0:N))
  allocate(a(N,N))
  allocate(atransp(N,N))

  x(0) = 3.0
  
  do i = 1, N
     x(i) = g(x(i-1))*(1-g(x(i-1)))-sqrt(g(x(i-1)))
     !  write(*,*) i, x(i)
  end do

  open(unit=123,file='saida.dat',status='unknown')
  do i=0,N,5
     write(123,*) i, x(i)
  end do
  close(unit=123)

  do i=1,N
     do j=1,N
        a(i,j) = sin(2.0*i)*cos(4.0*j*j)+1.0
     end do
  end do

  call transposta(N,a,atransp)

  do i=1,N
     do j=1,N
        write(*,*) i,j,a(i,j), atransp(i,j)
     end do
  end do
     
  deallocate(x)
  deallocate(a)
  deallocate(atransp)
  
contains

  function g(x)

    implicit none

    real :: x,g
    
    g = tan(x)-x*x
        
    return 
    
  end function g

  subroutine transposta(N,a,b)

    implicit none

    real, dimension(N,N) :: a, b
    integer :: i, j, N
    
    do i=1,N
       do j=1,N
          b(i,j) = a(j,i)
       end do
    end do
    
    return
    
  end subroutine transposta
  
end program vetores

