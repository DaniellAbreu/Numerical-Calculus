program primeiro

  implicit none

  integer :: i
  real :: x
  double precision :: y
  
  !do i=1,10

     !x = 1.0*i !real(i)
     
     !write(*,*) 'Olá número', i
     !write(*,*) 'O quadrado deste número é', i*i
     !write(*,*) 'A raiz quadrada deste número é', i**(0.5)
     !write(*,*) 'A raiz quadrada deste número é', sqrt(x)
     !write(*,*) 'Soma, subtração e divisão', i+1, x-0.5, i/x
     !write(*,*) 'divisão trivial?', 10/i, 10/x
     !write(*,*) 'precisão é importante:', 2.0d0 + 1.0d-9

     !write(*,*) mod(i,2)
     !if(mod(i,2) .eq. 0) then    ! .eq. , .ne. , .lt , .le. , .gt. , .ge.
     !   ==, /=, >=, <, <=...
     !   write(*,*) 'O número é par', i
     !else
     !   write(*,*) 'O número é ímpar', i
     !end if
 
  open(unit=123,file='saida.dat',status='unknown')
  
  do i = 0, 100

     write(123,*) i*2.0*acos(-1.0)/100.0, sin(i*2.0*acos(-1.0)/100.0), cos(i*2.0*acos(-1.0)/100.0)
        
  end do

  close(unit=123)
     
  !end do
  
end program primeiro

