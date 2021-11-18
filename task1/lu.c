#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../matrix.h"

int main()
{ // Solution by lu factorization algorithm
    printf("Solution of problem 4.37 with LU decomposition algorithm\n\n");

    matrix a = import_data("a.csv");
    printf("A = \n");
    print_matrix(a, 1);

    matrix b = import_data("b.csv");
    printf("B = \n");
    print_matrix(b, 1);

    clock_t tic = clock(); // Lets count how much time the calculation takes

    matrix *lup = lu_decompose(a);
    char *foo[] = {"L", "U", "P"};
    for (size_t i = 0; i < 3; i++)
    {
        printf("%s =\n", foo[i]);
        print_matrix(lup[i], 1);
        printf("\n");
    }

    /* The system is given by [a][x] = [b]
       where [x] is the solution of [u][x] = [y]
       and [y] is the solution of [l][y] = [p][b] */

    // Obtaining [y]
    matrix new_b = matrix_product(lup[2], b);
    matrix y = solve_lower(lup[0], new_b);

    // Obtaining [x]
    matrix x = solve_upper(lup[1], y);

    clock_t toc = clock();

    printf("Solution =\n");
    print_matrix(x, 12);
    printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

    free_matrix(&a);
    free_matrix(&b);
    free_matrix(&y);
    free_matrix(&x);

    for (size_t i = 0; i < 3; i++)
        free_matrix(&lup[i]);

    free(lup);

    return 0;
}