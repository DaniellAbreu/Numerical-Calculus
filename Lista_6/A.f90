program A

	implicit none
	
	double precision, dimension(:,:), allocatable :: Vander, VanderTransp, Gram
	double precision, dimension(:), allocatable :: t, y, aty, coefs, Ax, e
	integer :: i, j, k, Grau, N
	double precision :: Soma, erro_global
	
	! (Conjunto de N + 1 Pontos) e (Ajuste de Grau Determinado):
	N = 7
	Grau = 5
	
	allocate(Vander(N + 1, Grau + 1))
	allocate(VanderTransp(Grau + 1, N + 1))
	allocate(Gram(Grau + 1, Grau + 1))
	allocate(t(N + 1))
	allocate(y(N + 1))
	allocate(aty(Grau + 1))
	allocate(coefs(Grau + 1))
	allocate(Ax(N + 1))
	allocate(e(N + 1))
	
	! Dados Primeira Alternativa:
  open(unit=123,file='dados.dat', status='old', action='read')
		do i = 1, N + 1
			read(123,*) t(i), y(i)
		end do
	close(unit=123)
  
  ! Impressão dos Dados:
  do i = 1, N + 1
		write(*,*) t(i), y(i)
  end do
  
  ! Determinando Matriz de Vandermonde:
  do i = 1, N + 1
		do j = 1, Grau + 1
			Vander(i,j) = t(i) ** (j - 1)
		end do
  end do
  
  do i = 1, N + 1
		write(*,*) (Vander(i,j), j = 1, Grau + 1)
  end do
  
  ! Determinando a Matriz Transposta de Vandermonde:
  do i = 1, N + 1
		do j = 1, Grau + 1
			VanderTransp(j,i) = Vander(i, j)
		end do
  end do
  ! Possível Utilizar Função do Fortran: transpose().
  
  ! Determinando a Matriz de Gramm com a Multiplicação das Matrizes (Regularização do Sistema):
  do i = 1, Grau + 1
		do j = 1, Grau + 1
			Soma = 0.d0
			do k = 1, N + 1
				Soma = Soma + VanderTransp(i,k) * Vander(k,j)
			end do
			Gram(i,j) = Soma
		end do
  end do 
  
  ! Biblioteca OpenBLAS Otimizada: 
  ! call dgemm('N', 'N', Grau + 1, N + 1, 1.0d0, Vander, N + 1, VanderTransp, Grau + 1, 0.0d0, Gramm, Grau + 1)
	
	! Determinando a Multiplicação da Transposta com o Vetor Solução (Regularização do Sistema):
	do i = 1, Grau + 1
		Soma = 0.d0
		do j = 1, N + 1
			Soma = Soma + VanderTransp(i,j) * y(j)
		end do
		aty(i) = Soma
	end do
	
	! Solução do Sistema:
	call Gauss_Seidel(Gram, aty, coefs, Grau)
	
	! Cálculo Polinômio Interpolador:
	!call Gauss_Seidel(Vander, y, coefs, Grau)
	
	! Impressão de Coeficientes:
	write(*,*) coefs
	
	! Determinando a Multiplicação da Vandermonde Pelo Vetor de Coeficientes (||e|| = ||Ax - b||):
	do i = 1, N + 1
		Soma = 0.d0
		do j = 1, Grau + 1
			Soma = Soma + Vander(i,j) * coefs(j)
		end do
		Ax(i) = Soma
		e(i) = Ax(i) - y(i)
	end do
	
	! Calculo da Norma do Erro (||e||):
	erro_global = dsqrt(sum(e ** 2))
	
	! Imprimindo Erro Global:
	write(*,*) erro_global
	
	! Gerando Arquivo de Coeficientes:
	open(unit=123, file='coefs.dat', status='replace', action='write')
		write(123,'( *(G0, 1X) )') coefs
	close(unit=123)
	
	! Chamando Script de Geração de Gráfico:
	call system('gnuplot -persist plot_fit.gp')
	
	deallocate(Vander)
	deallocate(VanderTransp)
	deallocate(Gram)
	deallocate(t)
	deallocate(y)
	deallocate(aty)
	deallocate(coefs)
	deallocate(Ax)
	deallocate(e)
	
	contains
	
	subroutine Gauss_Seidel(A, b, x, Grau)
	
		implicit none
	
		double precision, dimension(Grau + 1, Grau + 1) :: A
		double precision, dimension(Grau + 1) :: b, x
		integer :: i, j, k, Grau
		double precision :: Soma
		
		do k = 1, 10000
			do i = 1, Grau + 1
				Soma = 0.d0
				do j = 1, Grau + 1
					if (j.ne.i) then
						Soma = Soma + x(j) * A(i,j)
					end if
				end do
				x(i) = (b(i) - Soma) / A(i,i)
			end do
		end do
		
		return
	
	end subroutine Gauss_Seidel
	
end program A
