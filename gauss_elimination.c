#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "matrix.h"

void gauss_eliminate(matrix *mtx, size_t args, ...)
{ /* Apply the Gauss Elimination algorithm to the matrix mtx and
    to all the vectors passed as arguments
    Note that all the matrixes' data will be changed
    If that is not intended, make sure to pass copies using copy_matrix() */

    va_list ap;
    int degenerate = 0;
    va_start(ap, args);

    for (size_t i = 0; i < args; i++)
    {                                     // Check if the number of rows of the arguments
        matrix *b = va_arg(ap, matrix *); // is equal to the number of rows of mtx
        if (b->rows != mtx->rows || b->cols != 1)
            error("The number of rows in the arguments of gauss_eliminate() must be equal");
    }

    va_end(ap);
    size_t h = 0;
    size_t k = 0;

    while (h < mtx->rows && k < mtx->cols)
    {

        if (fabs(mtx->elem[h][k]) < 0.05 && h != mtx->rows - 1)
        { // If pivot is null,
            size_t i_max;
            double max = 0;
            for (size_t n = h + 1; n < mtx->rows; n++)
            {
                if (fabs(mtx->elem[n][k]) >= max)
                {                                // evaluate the index of the element
                    max = fabs(mtx->elem[n][k]); // with the gratest absolute value in the rest
                    i_max = n;                   // of the column
                }
            }
            if (mtx->elem[i_max][k] == 0 || max <= fabs(mtx->elem[h][k]))
            { // If there is no new pivot, go to the next column
                degenerate = 1;
                k++;
                continue;
            }
            swap_rows(mtx, h, i_max); // If there is a pivot candidate, take it as the new pivot
            va_start(ap, args);
            for (size_t i = 0; i < args; i++)
            { // Swap the same rows in all the arguments
                matrix *b = va_arg(ap, matrix *);
                swap_rows(b, h, i_max);
            }
            va_end(ap);
        }

        for (size_t i = h + 1; i < mtx->rows; i++)
        { // Gauss elimination algorithm
            double m = mtx->elem[i][k] / mtx->elem[h][k];
            mtx->elem[i][k] = 0;
            for (size_t j = k + 1; j < mtx->cols; j++)
                mtx->elem[i][j] -= mtx->elem[h][j] * m;
            va_start(ap, args);
            for (size_t l = 0; l < args; l++)
            {
                matrix *b = va_arg(ap, matrix *);
                b->elem[i][0] -= b->elem[h][0] * m;
            }
            va_end(ap);
        }
        h++;
        k++;
    }

    if (degenerate)
        fprintf(stderr, "Warning: Near degenerate matrix detected. Results may not be precise.\n");
}

matrix solve(matrix *mtx, matrix *b)
{ // Return the solution x of a triangular system mtx * x = b
    if (mtx->rows != b->rows)
        error("The number of rows in the arguments of solve() must be equal");

    if (b->cols != 1)
        error("The second argument of solve() must be a vector");

    if (mtx->rows != mtx->cols)
        error("The matrix must be square and upper triangular in solve()");

    for (size_t j = 0; j < mtx->cols; j++)
    {
        for (size_t i = j + 1; i < mtx->rows; i++)
        {
            if (mtx->elem[i][j])
                error("The matrix must be upper triangular in solve()");
        }
    }

    size_t i = 0;
    while (i < mtx->rows)
    {
        if (null_row(*mtx, i))
        {
            if (b->elem[i][0])
            {
                printf("The system does not have a solution\n");
                return copy_matrix(*mtx);
            }
            else
            {
                remove_row(b, i);
                remove_row(mtx, i);
                continue;
            }
        }
        i++;
    }

    i = 0;
    while (i < mtx->cols)
    {
        if (null_col(*mtx, i))
        {
            remove_col(mtx, i);
            continue;
        }
        i++;
    }

    if (mtx->rows != mtx->cols)
    {
        printf("The system has infinitelly many solutions\n");
        return copy_matrix(*mtx);
    }

    int n = mtx->rows - 1;
    matrix sol = create_matrix(mtx->cols, 1);
    double sum;
    sol.elem[n][0] = b->elem[n][0] / mtx->elem[n][n];
    for (int k = n - 1; k >= 0; k--)
    {
        sum = 0;
        for (int j = k + 1; j <= n; j++)
            sum += mtx->elem[k][j] * sol.elem[j][0];
        sol.elem[k][0] = (b->elem[k][0] - sum) / mtx->elem[k][k];
    }
    return sol;
}

double det(matrix mtx)
{
    if (mtx.rows != mtx.cols)
        error("Argument is not a square matrix in det()");

    size_t h = 0;
    size_t k = 0;
    double det = 1;

    while (h < mtx.rows && k < mtx.cols)
    {

        if (mtx.elem[h][k] == 0 && h != mtx.rows - 1)
        {
            size_t i_max;
            double max = 0;
            for (size_t n = h + 1; n < mtx.rows; n++)
            {
                if (fabs(mtx.elem[n][k]) >= max)
                {
                    max = fabs(mtx.elem[n][k]);
                    i_max = n;
                }
            }
            if (mtx.elem[i_max][k] == 0)
            {
                k++;
                continue;
            }
            swap_rows(&mtx, h, i_max);
            det *= -1;
        }

        for (size_t i = h + 1; i < mtx.rows; i++)
        {
            double m = mtx.elem[i][k] / mtx.elem[h][k];
            mtx.elem[i][k] = 0;
            for (size_t j = k + 1; j < mtx.cols; j++)
                mtx.elem[i][j] -= mtx.elem[h][j] * m;
        }
        h++;
        k++;
    }

    for (size_t i = 0; i < mtx.rows; i++)
        det *= mtx.elem[i][i];

    return det;
}