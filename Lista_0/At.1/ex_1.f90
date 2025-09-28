program test

	implicit none
	
	integer :: i
	integer :: N
	real :: x
		
	open(unit=123, file="ex_1.dat", status="unknown")
	
	N = 25
	x = 0.0
	
	write(123,*) x, sin(x)
	
	do i = 1, N-1
		
		x = x + (2.0*acos(-1.0))/24
		write(123,*) x, sin(x)
		
	end do
		
	close(unit=123)

end program test
