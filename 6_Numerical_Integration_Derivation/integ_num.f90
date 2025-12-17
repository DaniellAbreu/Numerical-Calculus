program integral

  implicit none

  double precision :: h, soma, simp, retan
  integer :: N, i, k 

  open(unit=123,file='saida.dat',status='unknown')

  do k = 1,10
     N = 10*k 
     h = 2.d0/dble(N)

     !retângulos pela esquerda
     soma = 0.d0
     do i = 0,N-1
        soma = soma + y(i*h)*h
     end do
     retan = soma
     !write(*,*) 'Retângulo = ', soma


     !Simpson
     soma = 0.d0
     do i = 1,N/2
        soma = soma + y((2*i-1)*h)
     end do
     simp = 4.d0*soma
     soma = 0.d0
     do i = 1, (N-2)/2
        soma = soma + y(2*i*h)
     end do
     simp = simp + 2.d0*soma
     simp = (simp + y(0.d0)+y(N*h))*h/3.d0

     write(123,*) h, retan, simp
     
  end do

  

  close(unit=123)

   contains

     function y(t)

       implicit none

       double precision :: y, t

       !y = t**3+6.d0*t*t
       y = exp(-t*t)*sin(t)

       return

     end function y


   end program integral

