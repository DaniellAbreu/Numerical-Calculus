program diffin

  implicit none

  double precision, dimension (:), allocatable :: t, y, ylinha
  integer :: i, N, k
  double precision :: h, pi, ylinhaav, ylinhaat,ylinhace, ylinhaass

  open(unit=124,file='deriv.dat', status= 'unknown')
  
 ! do k = 0, 20
  k = 3
     N = 8*(2**k)
     pi = dacos(-1.d0)
     
     allocate(t(0:N))
     allocate(y(0:N))
     allocate(ylinha(0:N))
     
     h = 2.d0*pi/dble(N)
     do i = 0, N
        t(i) = i*h
        y(i) = dexp(-(t(i)**2)/23.d0)*datan(t(i)+0.1d0)*(dsqrt(t(i)))
     end do

     ylinha(0) = 0.5d0*(-3.d0*y(0)+4.d0*y(1)-y(2))/h
     do i=1,N-1
        ylinha(i) = 0.5d0*(y(i+1)-y(i-1))/h
     end do
     ylinha(N) = -0.5d0*(-3.d0*y(N)+4.d0*y(N-1)-y(N-2))/h
     
!!$     !Cálculo da derivada
!!$     i = N/8
!!$     ylinhaav = (y(i+1) - y(i))/h
!!$     ylinhaat = (y(i) - y(i-1))/h
!!$     ylinhace = 0.5d0*(y(i+1)-y(i-1))/h
!!$     ylinhaass = 0.5d0*(-3.d0*y(i)+4.d0*y(i+1)-y(i+2))/h
!!$     write(124,*) h, dabs(ylinhaav-dcos(pi/4.d0)), dabs(ylinhaat-dcos(pi/4.d0)), dabs(ylinhace-dcos(pi/4.d0)) &
!!$          , dabs(ylinhaass-dcos(pi/4.d0))
     !write(124,*) h, dabs(ylinhaav-dcos(0.d0)), dabs(ylinhaass-dcos(0.d0))
     
     open(unit=123,file='saida.dat',status='unknown')
     do i=0,N
        write(123,*) t(i), y(i), ylinha(i)
     end do
     close(unit=123)
     
     
     deallocate(t)
     deallocate(y)
     deallocate(ylinha)
 ! end do

  close(unit=124)
  
end program diffin


!Códigos gnuplot
!plot 'saida.dat' u 1:2 w p ps 1.5 pt 6, sin(x) w l 

