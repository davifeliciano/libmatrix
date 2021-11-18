#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../matrix.h"

int main()
{ // Solution by gauss elimination algorithm
    printf("Solution of problem 4.37 with Gauss Elimination algorithm\n\n");

    matrix a = import_data("a.csv");
    printf("A = \n");
    print_matrix(a, 1);

    matrix b = import_data("b.csv");
    printf("B = \n");
    print_matrix(b, 1);

    clock_t tic = clock(); // Lets count how much time the calculation takes

    gauss_eliminate(&a, 1, &b);

    printf("\nAfter the elimination, A and B will become\n");
    printf("A =\n");
    print_matrix(a, 1);

    printf("B =\n");
    print_matrix(b, 1);

    matrix sol = solve(&a, &b);

    clock_t toc = clock();

    printf("Solution by GE =\n");
    print_matrix(sol, 12);
    printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

    free_matrix(&a);
    free_matrix(&b);
    free_matrix(&sol);

    return 0;
}