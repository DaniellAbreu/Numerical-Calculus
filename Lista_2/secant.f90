program secant
	
	! Método da Secante - "Newton-Raphson" Sem Cálculo de Derivada - (Convergência Um Pouco Mais Lenta).
	
	implicit none
	
	double precision, dimension(:), allocatable :: x
	double precision :: tol
	integer :: N, M
	
	M = 100
	tol = 1.d-15
	allocate(x(-1:M))
	
	! Chutes Iniciais - (Necessário Dois Chutes Iniciais)
	x(-1) = 1.d0
	x(0) = 2.d0
	
	N = 0
	
	do while (N < size(x) - 2)
		x(N + 1) = x(N) - f(x(N)) * ((x(N) - x(N - 1)) / (f(x(N)) - f(x(N - 1))))
		
		if (abs(x(N + 1) - x(N)) < tol) then
			write(*,*) 'Convergiu para: ', x(N + 1)
			stop
		end if
		
		write(*,*) N, x(N + 1)
		N = N + 1
	end do
	
	deallocate(x)
	
	contains
	
	function f(x)
	
	implicit none
	
	double precision :: f, x
	
	f = (x ** 3.d0) - x - 1.d0
	
	return
	
	end function
	
end program secant
