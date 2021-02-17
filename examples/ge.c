#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../matrix.h"

int main(int argc, char **argv)
{ /* This program solver a linear system [a][x] = [b],
    [a] and [b] are imported from two csv files given as cl parameters */
    if (argc != 3)
    {
        fprintf(stderr, "Usage example: ./ge a.csv b.csv\n");
        exit(EXIT_FAILURE);
    }

    matrix a = import_data(argv[1]);
    printf("A = \n");
    print_matrix(a, 1);

    matrix b = import_data(argv[2]);
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