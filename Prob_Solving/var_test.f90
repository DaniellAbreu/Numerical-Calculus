program var_test
    implicit none
    integer :: x,y

    print *, 'Enter two natural numbers: '
    read (*,*) x,y

    print *, 'The sum and the multiplication of x and y are: ', x+y, x*y

end program var_test