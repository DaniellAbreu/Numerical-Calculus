program A
  
  ! Solução de Sistemas Algébricos Lineares - Processos Iterativos Matriciais.
  ! Eliminação Gaussiana; Pivoteamento; Substituição Reversa e Rotina de Verificação.
  
  implicit none
	
  integer :: i, j, k, N
  double precision, dimension(:,:), allocatable :: M
  double precision, dimension(:), allocatable :: Var, Temp
  double precision :: Lambda, Soma
	
  N = 3 ! Grau do Sistema Linear.
	
  allocate(M(N,N + 1)) ! Matriz Aumentada (Matriz Quadrada + Solução).
  allocate(Var(N))

  M(1,1) = 2.d0
  M(1,2) = -1.5d0
  M(1,3) = 3.d0
  M(1,4) = 1.d0
  
  M(2,1) = -1.d0
  M(2,2) = 0.d0
  M(2,3) = 2.d0
  M(2,4) = 3.d0
  
  M(3,1) = 4.d0
  M(3,2) = -4.5d0
  M(3,3) = 5.d0
  M(3,4) = 1.d0
	
	write(*,*) 'Antes do Escalonamento'
  do i = 1, N
    write(*,*) (M(i,j),j = 1, N + 1)
  end do
	
	do j = 1, N ! Percorre Colunas.
		if (M(j,j) == 0) then ! Método do Pivoteamento.
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
					write(*,*) 'Problema! A Matriz é Singular!' ! Determinante(M) = 0 / Não é Invertível.
					stop
				end if
			end do
		else
			do i = j + 1, N ! Percorre Linhas Abaixo dos Pivôs.
				Lambda = M(i,j)/M(j,j) ! Relação de Substituição de Linhas (Elemento Atual / Pivô).
				
				do k = 1, N + 1
					M(i,k) = M(i,k) - Lambda * M(j,k) ! Troca dos Elementos da Linha (Cuidar Com Matriz Aumentada).
				end do
			end do
		end if
	end do

	write(*,*) 'Depois do Escalonamento - Upper Triangular'
  do i = 1, N
    write(*,*) (M(i,j),j = 1, N + 1)
  end do
  
  ! Algorítmo de Substituição Reversa.
  Var(N) = M(N,N + 1)/M(N,N) ! Solução do Última Variável.
  do i = N - 1, 1, -1 ! Percorre Linhas Decrescentemente.
		Soma = 0.d0
		do k = i + 1, N
			Soma = Soma + M(i,k) * Var(k) ! Soma das Variáveis Resolvidas.
		end do
		Var(i) = (M(i, N + 1) - Soma)/M(i,i) ! Solução da Variável Atual.
  end do
  
  write(*,*) 'Solução do Sistema Linear'
  do i = 1, N
		write(*,*) i, Var(i)
  end do
  
	! Rotina de Verificação.
	Soma = 0.d0
	do i = 1, N
		Soma = Soma + M(1,i) * Var(i)
	end do
	
	if (Soma == M(1, N + 1)) then
		write(*,*) 'Solução Verificada com Sucesso! :)'
	else
		write(*,*) 'Solução Incorreta :(', Soma, M(1, N + 1)
	end if

	deallocate(M)
	deallocate(Var)
	
end program A
