program float
  use, intrinsic :: iso_c_binding, only: sp=>c_float, dp=>c_double
  implicit none

  real(sp) :: float32
  real(dp) :: float64

  float32 = 1.45
  float64 = 43.65

  print*, float32, float64

end program float