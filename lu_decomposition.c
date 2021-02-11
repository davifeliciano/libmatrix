#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

matrix *lu_decompose(matrix mtx)
{ // Return a pointer to an array with the matrixes l, u and p
    if (mtx.rows != mtx.cols)
        error("The matrix must be square in lu_decompose()");

    for (size_t i = 0; i < mtx.rows; i++)
    {
        if (null_row(mtx, i))
            error("There is a null row in the argument of lu_decompose()");
        if (null_col(mtx, i))
            error("There is a null column in the argument of lu_decompose()");
    }

    matrix l = create_matrix(mtx.rows, mtx.cols);
    matrix u = create_matrix(mtx.rows, mtx.cols);
    matrix p = Id(mtx.rows); // Permutation matrix, for pivoting
    print_matrix(p, 1);

    int degenerate = 0;

    for (size_t i = 0; i < mtx.rows; i++)
    {
        if (fabs(mtx.elem[i][i] < 0.05))
        { // If pivot is null,
            size_t i_max;
            double max = 0;
            for (size_t n = i + 1; n < mtx.rows; n++)
            {
                if (fabs(mtx.elem[n][i]) >= max)
                {                               // evaluate the index of the element
                    max = fabs(mtx.elem[n][i]); // with the gratest absolute value in the rest
                    i_max = n;                  // of the column
                }
            }
            if (mtx.elem[i_max][i] && max > mtx.elem[i][i])
            { // If there is a new pivot, swap rows in mtx and in the permutation matrix p
                swap_rows(&mtx, i, i_max);
                swap_rows(&p, i, i_max);
            }
            else
                degenerate = 1; // If no new pivot, matrix is near degenerate
        }
    }

    if (degenerate)
        fprintf(stderr, "Warning: Near degenerate matrix detected. Results may not be precise.\n");

    // lu decomposition algorithm
    for (size_t i = 0; i < mtx.rows; i++)
    {
        l.elem[i][0] = mtx.elem[i][0];
        u.elem[i][i] = 1;
    }

    for (size_t j = 1; j < mtx.cols; j++)
        u.elem[0][j] = mtx.elem[0][j] / l.elem[0][0];

    for (size_t i = 1; i < mtx.rows; i++)
    {
        for (size_t j = 1; j <= i; j++)
        {
            double sum = 0;
            for (size_t k = 0; k < j; k++)
                sum += l.elem[i][k] * u.elem[k][j];

            l.elem[i][j] = mtx.elem[i][j] - sum;
        }

        for (size_t j = i + 1; j < mtx.rows; j++)
        {
            double sum = 0;
            for (size_t k = 0; k < i; k++)
                sum += l.elem[i][k] * u.elem[k][j];

            u.elem[i][j] = (mtx.elem[i][j] - sum) / l.elem[i][i];
        }
    }

    matrix *result = (matrix *)malloc(3 * sizeof(matrix));
    result[0] = l;
    result[1] = u;
    result[2] = p;
    return result;
}