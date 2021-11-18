#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <time.h>
#include "../matrix.h"

int main(int argc, char **argv)
{ // Decompose a matrix in Q R product
    if (argc != 2)
    {
        fprintf(stderr, "Usage: ./%s [filename]\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    matrix data = import_data(argv[1]);

    printf("Loaded:\n");
    print_matrix(data, 12);

    clock_t tic = clock(); // Lets count how much time the calculation takes
    matrix *qr = qr_decompose(data);
    clock_t toc = clock();

    matrix prod = matrix_product(qr[0], qr[1]);
    char *foo[] = {"Q", "R"};

    for (size_t i = 0; i < 2; i++)
    {
        printf("%s =\n", foo[i]);
        print_matrix(qr[i], 12);
        free_matrix(&qr[i]);
    }

    printf("Q * R =\n");
    print_matrix(prod, 12);
    printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

    free(qr);
    free_matrix(&data);
    free_matrix(&prod);
    return 0;
}