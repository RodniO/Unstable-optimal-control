# Unstable-optimal-control

Finds a curve of minimum length from (0;1) to (1;0) with integral of x*y under the curve equal to 0 (see the figure). Algorithm from
[1] Осинский А. И. ПРИБЛИЖЕННОЕ РЕШЕНИЕ ЧАСТНОЙ НЕУСТОЙЧИВОЙ ЗАДАЧИ ОПТИМАЛЬНОГО УПРАВЛЕНИЯ // International Conference "Optimal Control and Differential Games" dedicated to the 110th anniversary of L. S. Pontryagin, 2018. DOI: 10.4213/proc23022

![alt text](https://github.com/RodniO/[reponame]/blob/master/lsp.png?raw=true)

To compile and run it, use either

make gnu_run

(Requires gfortran, BLAS and LAPACK)

or

make run

(Requires ifort and mkl)