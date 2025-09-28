program ex_4

	implicit none
	
	real :: r
	real :: i
    call random_number(r)
    
    i = (r * 101)
    write(*,*) i

end program ex_4
