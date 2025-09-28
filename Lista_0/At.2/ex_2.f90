program ex_2

	implicit none
	
	integer :: i
	integer, parameter :: N = 25
	real :: data(N, 2)
	real :: new_data(N)
		
	open(unit=123, file="ex_1_copy.dat", status="old")
	
	do i = 1, N
		
		read(123, *) data(i, 1), data(i, 2)
		
	end do
		
	close(unit=123)
	
	do i = 1, N
		
		new_data(i) = asin(data(i, 2))
		
	end do
	
	open(unit=123, file="ex_1_copy.dat", status="replace")
	
	do i = 1, N
	
		write(123, *) data(i, 1), data(i, 2), new_data(i)
	
	end do
	
	close(unit=123)

end program ex_2
