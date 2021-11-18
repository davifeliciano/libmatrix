#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../matrix.h"

int main(int argc, char **argv)
{ // Compute the eigenvalues of a matrix with qr factorization method
    if (argc != 2)
    {
        fprintf(stderr, "Usage: ./%s [filename]\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    matrix data = import_data(argv[1]);

    printf("Loaded:\n");
    print_matrix(data, 12);

    clock_t tic = clock(); // Lets count how much time the calculation takes
    matrix eigen = eigenvalues_qr(data);
    clock_t toc = clock();

    printf("Real Eigenvalues:\n");
    print_matrix(eigen, 12);

    printf("\nElapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

    free_matrix(&data);
    free_matrix(&eigen);

    return 0;
}