This project evaluates numerically the solution of the problem 4.37 of the book
Numerical Methods for Engineers and Scientists: An Introduction with Applications using MATLAB
by Amos Gilat and Vish Subramaniam. The data of the linear system to be solved are stored in
two csv files, a.csv and b.csv. Take a look on the problem and check the values.

task1.nb is a Mathematica Notebook File that show the solution, for comparation purposes.
task1.png is a screenshot of this Mathematica Notebook.

The solution was obtained using two methods. Gaussian Elimination, in ge.c and LU Factorization
in lu.c. To compile both programs, type

    make all && make clean

in a linux terminal. Some examples of usage of this matrix library are shown in examples directory.
examples/ge.c and examples/lu.c are versions of these programs that take arbitrary matrixes as
csv files passed as command line arguments. To compile them type the same command as before, but in examples dir.
I have not intensivelly tested any of the algorithms, so keep in mind that you may find some bugs.

Author: Davi Feliciano      Email: dfeliciano37@gmail.com