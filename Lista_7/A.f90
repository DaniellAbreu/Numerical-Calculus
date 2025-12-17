program A

	implicit none
	
	double precision, dimension(:), allocatable :: x, y, ylinha_cent, ylinhalinha_cent
	double precision :: Ain, Bout, h
	integer :: i, Np
	
	! Intervalo da Função
	Ain = 0.d0
	Bout = 2 * dacos(-1.d0)
	
	! Número de Pontos
	Np = 21
	
	! Divisão do Intervalo da Função em (Np - 1) Sub-Intervalos.
	h = dabs(Bout - Ain) / dble(Np - 1) ! Intervalo de Discretização.
	
	allocate(x(0:Np - 1))
	allocate(y(0:Np - 1))
	allocate(ylinha_cent(Np - 2))
	allocate(ylinhalinha_cent(Np - 2))
	
	do i = 0, Np - 1
		x(i) = i * h
	end do
	
	do i = 0, Np - 1
		y(i) = f(x(i))
	end do
	
	do i = 1, Np - 2
		! Cálculo por Diferenças Finitas Centradas da Primeira Derivada da Função (Segunda Ordem):
		ylinha_cent(i) = 0.5d0 * (f(x(i + 1)) - f(x(i - 1))) / h
		
		! Cálculo por Diferenças Finitas Centradas da Segunda Derivada da Função (Segunda Ordem):
		ylinhalinha_cent(i) = (f(x(i + 1)) - 2.d0 * f(x(i)) + f(x(i - 1))) / h**2
	end do 

	open(unit=123, file="sine_function.dat", status="unknown")
		do i = 0, Np - 1
			write(123,*) x(i), y(i)
		end do 
	close(123)
	
	open(unit=124, file="derivatives.dat", status="unknown")
		do i = 1, Np - 2
			write(124,*) x(i), ylinha_cent(i), dcos(x(i)), ylinhalinha_cent(i), -dsin(x(i))
		end do 
	close(124)
	
	deallocate(x)
	deallocate(y)
	deallocate(ylinha_cent)
	deallocate(ylinhalinha_cent)
	
	contains
	
	function f(x)

	implicit none
	
	double precision :: f, x
	
	f = dsin(x)
	
	return
	
	end function
	
end program A
