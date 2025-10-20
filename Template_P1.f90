program sist_naolin

  implicit none

  double precision, dimension (:,:), allocatable :: a
  double precision, dimension (:), allocatable :: x, y, b
  integer :: n, Nfim, size
  
  size = 2 !número de equações'
  Nfim = 500

  !alocação das variáveis no tamanho desejado
  allocate(x(size))
  allocate(y(size))
  allocate(b(size))
  allocate(a(size,size))

  !preenchimento aleatório das variáveis vetoriais/matriciais, apenas para fins de compilação -- alterar para a necessidade específica do programa.
  a = 1.d0
  b = 1.d0
  x = 1.d0
  y = 0.d0

  !Laço para realizar um processo iterativo com a função f
  do n = 1, Nfim
     x = f(size,x)
     write(*,*) n, x
  end do

  !Laço para chamar o algoritmo de solução de sistema por Nfim vezes utilizando uma subrotina
  do n = 1, Nfim
     call resolve_sis(size,a,f(size,x),y) !exemplo de chamamento da subrotina envolvendo o uso de funções. A ser alterado de acordo com o uso a ser dado à subrotina.
  end do

  deallocate(x)
  deallocate(y)
  deallocate(b)
  deallocate(a)
    
contains

  function f(size,x) !Definição de uma função vetorial. Para definir uma função escalar, remover a variável "size".

    implicit none

    double precision, dimension (size) :: x, f
    integer :: size

    !se a função for escalar, remover tudo abaixo e escrever apenas
    !f = dsin(x) + 1.d0 + dcos(x) !exemplo aleatório de função escalar.
    
    f(1) = 1.d0+dcos(x(1))+dsin(x(2)) !função aleatória apenas para fins compilação. Substituir pela função relevante. Neste caso, tomou-se size = 2.
    f(2) = dlog(dabs(x(2))+1.d0)

    return

  end function f

  subroutine resolve_sis(size,a,b,x)

    implicit none

    double precision, dimension (size,size) :: a !matriz dos coeficienrtes do sistema
    double precision, dimension (size) :: b,x !vetor solução e vetor resposta do sistema
    integer :: size !tamanho do sistema
    integer :: i, j, k  !variáveis auxiliares para implementar o algoritmo
    double precision :: soma !variáveis auxiliares para implementar o algoritmo
    
    do k = 1, 3000 !repetir por 3000 vezes o processo iterativo
       x = f(size,b) !função arbitrária, apenas para efeitos de compilação do programa -- linha a ser removida antes do uso do programa.
       !incluir o algortimo do processo iterativo aqui. Remover a linha anterior
       !a variável de saída desta subrotina é a variável x. Atenção à orde de preenchimento das variáveis!
    end do

    return

  end subroutine resolve_sis


end program sist_naolin

