program Matr_iter
  
  ! Processos Iterativos Matriciais Lineares - Assumindo Matrizes Quadradas.
  ! Método de Gauss-Seidel, Gauss-Jacobi, Resíduo e SOR.
  
  ! !!!IMPORTANTE!!!
  
  ! Processo Iterativo: Ax = b | x(n + 1) = Tx(n) + c.
  
  ! Matriz 'A' é Invertível -> Possui Solução Única.
  
  ! Ponto Fixo é Assintoticamente Estável <=> O Módulo de TODOS os Autovalores da Matriz 'T' são < 1 | (Raio Espectral < 1).
  ! Por Consequência, se Satisfeito a Condição Acima <=> A Matriz 'T' é Convergente.
  ! Polinômio Característico: (T - Lambda * I) * v = 0.
  
  ! Raio Espectral = Autovalor de Maior Módulo da Matriz 'T'.
  ! Raio Espectral Dita o Quão Rápido o Processo Iterativo irá Convergir - (Menor Raio Espectral -> Menos Iterações).
  
  ! Matriz Estritamente Diagonalmente Dominante é Convergente Para Ambos os Métodos | (Raio Espectral < 1).
  
  implicit none

  integer :: i, j, k, N
  double precision :: soma, tol, x_sei, x_old, Omega ! Omega = Parâmetro de Relaxação.
  double precision, dimension (:,:), allocatable :: A
  double precision, dimension (:), allocatable :: b, x, xold

  N = 5
  tol = 1.d-13
  
  ! Sub-Relaxado | 0 < Omega < 1.
  ! Gauss-Siedel | Omega = 1.
  ! Sobre-Relaxado | 1 < Omega < 2.
  Omega = 1.d0
  
  allocate(A(N,N))
  allocate(b(N))
  allocate(x(N))
  allocate(xold(N))

  A(1,1) = 8.d0
  A(1,2) = 1.d0
  A(1,3) = -1.d0
  A(1,4) = 1.d0
  A(1,5) = 3.d0
  b(1) = 7.d0

	A(2,1) = 3.d0
  A(2,2) = 18.d0
  A(2,3) = -4.d0
  A(2,4) = 0.d0
  A(2,5) = 5.d0
  b(2) = 6.d0
  
  A(3,1) = 0.d0
  A(3,2) = -2.d0
  A(3,3) = 7.d0
  A(3,4) = 1.d0
  A(3,5) = -1.d0
  b(3) =  2.d0

  A(4,1) = 1.d0
  A(4,2) = 0.d0
  A(4,3) = 2.d0
  A(4,4) = -8.d0
  A(4,5) = 1.d0
  b(4) = 2.d0

  A(5,1) = 1.d0 
  A(5,2) = -1.d0 
  A(5,3) = -1.d0
  A(5,4) = -1.d0
  A(5,5) = 9.d0
  b(5) = 3.d0

  ! Método é Insensível ao Chute inicial | Não Afeta Convergência em Casos de Matriz Bem Comportada.
  ! Caracterização do Ponto Fixo é Bem Definida em Processos Lineares = Ax + b.
  
  ! Chutes Iniciais de Solução do Sistema.
  x(1) = 0.d0
  x(2) = 0.d0
  x(3) = 0.d0
  x(4) = 0.d0
  x(5) = 0.d0
  
  
  ! Notação: x = (x(1), x(2), x(3)) = (x, y, z).
	! Processo Iterativo de Gauss-Jacobi = x(n + 1) = [D ** (-1)] * [b - (L + U) * x(n)].
	! Processo Iterativo de Gauss-Siedel = x(n + 1) = [(D + L) ** (-1)] * [b - U * x(n)].
	! Gauss-Siedel = Atualiza a Variável a Partir da Melhor Estimativa Encontrada.
	
  k = 0
  do while (residuo(b, A, x, N) .ge. tol) ! Processo Iterativo.
     k = k + 1
     ! xold = x ! Cópia de x para Gauss-Jacobi.
     do i = 1, N ! Percorrer Todas as Linhas do Sistema.
        soma = 0.d0
        do j = 1, N ! Percorrer Todas as Colunas do Sistema.
           if(j.ne.i) then ! Apenas Termos Não-Diagonais.
              soma = soma + A(i,j) * x(j) ! xold para Gauss-Jacobi, x para Gauss-Seidel.
           end if
        end do
        x_old = x(i)
        x_sei = (b(i) - soma) / A(i,i)
        x(i) = (1 - Omega) * x_old + Omega * x_sei ! Processo Iterativo Por SOR - (Successive Over-Relaxation).
     end do
     write(*,*) k, residuo(b, A, x, N)
  end do
  write(*,*) 'A Solução é', x

  deallocate(A)
  deallocate(b)
  deallocate(x)
  deallocate(xold)

contains

  function residuo(b, A, x, N)

    implicit none

    integer :: i, N, j
    double precision, dimension (N,N) :: A
    double precision, dimension (N) :: b,x,c
    double precision :: residuo, soma
    
    do i = 1, N
       soma = 0.d0
       do j = 1, N
          soma = soma + A(i,j) * x(j)
       end do
       c(i) = b(i) - soma
    end do

    residuo = 0.d0
    do i = 1, N
       residuo = residuo + c(i) * c(i)
    end do
    residuo = dsqrt(residuo)
    
    return
    
  end function residuo
  
end program Matr_iter

