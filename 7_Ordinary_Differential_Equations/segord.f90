program segord

  implicit none

  double precision, dimension (:), allocatable :: y
  double precision :: h
  integer :: ord, Nfim, n 
  
  ord = 2

  Nfim = 64000
  h = 4.d0/dble(Nfim)
  
  allocate(y(ord))

  open(unit=123,file='saida.dat', status= 'unknown')
  
  y(1) = 1.d0
  y(2) = 0.d0
  write(123,*) 0.d0, y
  do n = 0, Nfim
     y = y + h*f(n*h,y,ord)
     write(123,*) (n+1)*h, y
  end do

  close(unit=123)
  
  deallocate(y)
  
contains

  function f(t,y,ord)

    double precision, dimension (ord) :: f, y
    double precision :: t
    integer :: ord
        
    f(1) = y(2)
    f(2) = -4.d0*y(1)

    return
    
  end function f
  

end program segord

