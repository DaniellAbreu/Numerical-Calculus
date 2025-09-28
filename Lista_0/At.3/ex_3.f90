program ex_3
    
    implicit none

    integer :: i
    integer, parameter :: N = 25
    real :: data(N, 2)
    real :: ref

    ! Read numbers from file
    open(unit=123, file="ex_1_copy.dat", status="old")
    
    do i = 1, N
        
        read(123, *) data(i, 1), data(i, 2)
    
    end do
    
    close(unit=123)
 
    open(unit=124, file="ex_3.dat", status="replace")
    
    call bubble_sort(data, N)
    
    do i = 1, N
    
        write(124, *) i, data(i, 2)
        
    end do
    
    close(unit=124)

	contains
	
	subroutine bubble_sort(data, N)
	
	implicit none
	
	integer :: i, S
	integer :: N
    real :: data(N, 2)
    real :: ref
	
	S = 1
    do while (S /= 0)
        S = 0
        do i = 1, N - 1
            if (data(i, 2) > data(i + 1, 2)) then
                ref = data(i, 2)
                data(i, 2) = data(i + 1, 2)
                data(i + 1, 2) = ref
                S = S + 1
            end if
        end do
    end do
	
	end subroutine bubble_sort
	
end program ex_3
