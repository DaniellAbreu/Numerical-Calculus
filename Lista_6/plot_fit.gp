reset
set title "Ajustes Polinomiais em Conjuntos de Dados"
set xlabel "x"
set ylabel "y"
set grid
set key left top

# Leitura de Coeficientes:

coefs = system("cat coefs.dat")
n = words(coefs)

# Ajuste Polinomial:

f(x) = sum [i=0:n-1] ( real(word(coefs, i+1)) * x**i )

plot f(x) with lines lw 2 title "Ajuste de Curva", "dados.dat" using 1:2 with points pt 7 title "Dados"
