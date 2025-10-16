program B
  
  ! Solução de Sistemas Algébricos Lineares - Processos Iterativos Matriciais.
  ! Eliminação Gaussiana; Pivoteamento Simples; Substituição Reversa e Rotina de Verificação.
  ! Fatoração LU.
  
  implicit none
	
  integer :: i, j, k, N
  double precision, dimension(:,:), allocatable :: M, U, Lambda
  double precision, dimension(:), allocatable :: C, X, Temp, D
  double precision :: Soma
	
  N = 4 ! Grau do Sistema Linear.
	
  allocate(M(N,N + 1)) ! Matriz Aumentada (Matriz Quadrada + Solução).
  allocate(X(N)) ! Valores de Variáveis Solução de U.
  allocate(C(N)) ! Valores de Variáveis Solução de L.
  allocate(D(N)) ! Vetor Solução do Sistema.
  
  ! Alocação de Matrizes para Fatoração LU.
  allocate(U(N,N))
  allocate(Lambda(N,N))

  M(1,1) = acos(-1.d0)
  M(1,2) = exp(1.d0)
  M(1,3) = sqrt(2.d0)
  M(1,4) = -sqrt(3.d0)
  M(1,5) = 1.d0
  
  M(2,1) = acos(-1.d0) ** 2
  M(2,2) = exp(1.d0)
  M(2,3) = -(exp(1.d0) ** 2)
  M(2,4) = 3.d0 / 7.d0
  M(2,5) = 1.d0
  
  M(3,1) = sqrt(5.d0)
  M(3,2) = -sqrt(6.d0)
  M(3,3) = 1.d0
  M(3,4) = -sqrt(2.d0)
  M(3,5) = 1.d0
  
  M(4,1) = acos(-1.d0) ** 3
  M(4,2) = exp(1.d0) ** 2
  M(4,3) = -sqrt(7.d0)
  M(4,4) = 1.d0 / 9.d0
  M(4,5) = 1.d0
  
  ! Formato Inicial da Matriz L - Fatoração LU.
  Lambda = 0.d0
  do i = 1, N
		Lambda(i,i) = 1.d0
  end do
  
  ! Alocação do Vetor Solução.
  do i = 1, N
		D(i) = M(i, N + 1)
  end do
	
	write(*,*) 'Antes do Escalonamento - Matriz Base'
  do i = 1, N
    write(*,*) (M(i,j), j = 1, N + 1)
  end do
	
	do j = 1, N ! Percorre Colunas.
		if (M(j,j) == 0) then ! Método do Pivoteamento Simples.
			do k = j + 1, N
				if (M(k,j) /= 0) then
					allocate(Temp(N + 1))
					do i = 1, N + 1
						Temp(i) = M(j,i)
						M(j,i) = M(k,i)
						M(k,i) = Temp(i)
					end do
					deallocate(Temp)
				else
					write(*,*) 'Problema! A Matriz é Singular!' ! Determinante(M) = 0 - Não é Invertível.
					stop
				end if
			end do
		else
			do i = j + 1, N ! Percorre Linhas Abaixo dos Pivôs.
				Lambda(i,j) = M(i,j)/M(j,j) ! Relação de Substituição de Linhas (Elemento Atual / Pivô).
				
				do k = 1, N + 1
					M(i,k) = M(i,k) - Lambda(i,j) * M(j,k) ! Troca dos Elementos da Linha (Cuidar Com Matriz Aumentada).
				end do
			end do
		end if
	end do

	! Criação da Matriz U - Fatoração LU.
	do i = 1, N
		do j = 1, N
			U(i,j) = M(i,j)
		end do
	end do

	write(*,*) 'Depois do Escalonamento: Matriz Base'
  do i = 1, N
    write(*,*) (M(i,j), j = 1, N + 1)
  end do
  
  write(*,*) 'Matriz U - Upper Triangular'
  do i = 1, N
    write(*,*) (U(i,j), j = 1, N)
  end do
  
  write(*,*) 'Matriz L - Lower Triangular'
  do i = 1, N
    write(*,*) (Lambda(i,j), j = 1, N)
  end do
	
	! Algoritmo de Substituição Direta
	C(1) = D(1)
	do i = 2, N
		Soma = 0.d0
		do k = 1, i - 1
			Soma = Soma + Lambda(i, k) * C(k)
			end do
		C(i) = D(i) - Soma
	end do
	
	! Algorítmo de Substituição Reversa.
  X(N) = C(N)/U(N,N) ! Solução do Última Variável.
  do i = N - 1, 1, -1 ! Percorre Linhas Decrescentemente.
		Soma = 0.d0
		do k = i + 1, N
			Soma = Soma + U(i,k) * X(k) ! Soma das Variáveis Resolvidas.
		end do
		X(i) = (C(i) - Soma) / U(i,i) ! Solução da Variável Atual.
  end do
  
  write(*,*) 'Solução do Sistema Linear'
  do i = 1, N
		write(*,*) i, X(i)
  end do
	
	! Rotina de Verificação.
	Soma = 0.d0
	do i = 1, N
		Soma = Soma + M(1,i) * X(i)
	end do
	
	if (Soma == M(1, N + 1)) then
		write(*,*) 'Solução Verificada com Sucesso! :)'
	else
		write(*,*) 'Solução Incorreta :(', Soma, M(1, N + 1)
	end if
	
	deallocate(M)
	deallocate(X)
	deallocate(C)
	deallocate(D)
	deallocate(U)
	deallocate(Lambda)
	
end program B
