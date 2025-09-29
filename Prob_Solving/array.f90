program array_slice
  implicit none

  integer :: i
  integer :: array1(10)  ! 1D integer array of 10 elements

  integer :: array2(10, 10)  ! 2D integer array of 100 elements

  array1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]  ! Array constructor
  print *, "Array at odd indices"
  print *, array1(1:10:2)  ! Print out elements at odd indices, start at element 1, go up to element 10, in strides of 2

  print *, "Array in reverse!"
  print *, array1(10:1:-1)  ! Print an array in reverse

  array1(:) = 0  ! Set all elements to zero
  print *, "Reset array to 0!"
  print *, array1

  array1 = [(i, i = 1, 10)]  ! Implied do loop constructor
  print *, "Now back to 1 to 10!"
  print *, array1


  array1(1:5) = 1  ! Set first five elements to one
  print *, "Set first five elements to 1"
  print *, array1

  array1(6:) = 2  ! Set all elements after five to two
  print *, "Set all elements after 5 to 2"
  print *, array1


  print *, "First column of 2d array"
  print *, array2(:,1)  ! Print out the first column in a 2D array


end program array_slice