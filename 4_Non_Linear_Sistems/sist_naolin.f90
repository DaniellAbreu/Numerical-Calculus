program sist_naolin

  implicit none

  double precision, dimension (:), allocatable :: x, y
  integer :: n, Nfim, size
  
  size = 2 !número de equações
  Nfim = 10
  
  allocate(x(size))
  allocate(y(size))

  x(1) = 7.d0
  x(2) = 3.d0
  
  do n = 1, Nfim
     call resolve_sis(size,gradf(size,x),f(size,x),y)
     x = x - y
     write(*,*) n, x
     !read(*,*)
  end do

  deallocate(x)
  deallocate(y)
  
contains

  function f(size,x)

    implicit none

    double precision, dimension (size) :: x, f
    integer :: size

    f(2) = x(1)*x(1)+x(2)*x(2)-64.d0
    f(1) = x(1)*x(2)-25.d0

    return

  end function f

  function gradf(size,x)

    implicit none

    double precision, dimension (size,size) :: gradf
    double precision, dimension (size) :: x  
    integer :: size

    gradf(2,1) = 2.d0*x(1)
    gradf(2,2) = 2.d0*x(2)
    gradf(1,1) = x(2)
    gradf(1,2) = x(1)

    return

  end function gradf

  subroutine resolve_sis(size,a,b,x)

    implicit none

    double precision, dimension (size,size) :: a
    double precision, dimension (size) :: b,x
    integer :: size, i, j, k
    double precision :: soma
    
    do k = 1, 3000
       do i = 1, size
          soma = 0.d0
          do j = 1, size
             if(j.ne.i) then
                soma = soma + a(i,j)*x(j)
             end if
          end do
          x(i) = (b(i)-soma)/a(i,i)
       end do
       !write(*,*) k, x
    end do

    return

  end subroutine resolve_sis


end program sist_naolin

