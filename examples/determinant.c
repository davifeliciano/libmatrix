#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../matrix.h"

int main(int argc, char **argv)
{ /* Calculate the determinant of a matrix given by csv file and time spent computing 
    Filename is given as the first cl argument */
    if (argc > 2)
    {
        fprintf(stderr, "Usage: %s [FILE]\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    matrix a = import_data(argv[1]);
    print_matrix(a, 5);
    printf("\n");

    clock_t tic = clock();
    printf("Determinant = %.3f\n", det(a));
    clock_t toc = clock();

    free_matrix(&a);
    fprintf(stderr, "Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
    return 0;
}