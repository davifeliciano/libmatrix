#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

matrix inverse(matrix mtx)
{ // Returns the inverse of mtx, without changing its elements
    if (mtx.rows != mtx.cols)
        error("Argument is not a square matrix in inverse()");

    matrix copy = copy_matrix(mtx);
    matrix inv = Id(copy.rows);
    size_t h = 0;
    size_t k = 0;

    while (h < copy.rows && k < copy.cols)
    {

        if (fabs(copy.elem[h][k]) < 1.0 && h != copy.rows - 1)
        {
            size_t i_max;
            double max = 0;
            for (size_t n = h + 1; n < copy.rows; n++)
            {
                if (fabs(copy.elem[n][k]) >= max)
                {
                    max = fabs(copy.elem[n][k]);
                    i_max = n;
                }
            }
            if (copy.elem[i_max][k] == 0)
            {
                k++;
                continue;
            }
            swap_rows(&copy, h, i_max);
            swap_rows(&inv, h, i_max);
        }

        double n = 1 / copy.elem[h][k];
        for (size_t j = 0; j < copy.cols; j++)
        {
            copy.elem[h][j] *= n;
            inv.elem[h][j] *= n;
        }

        for (size_t i = 0; i < copy.rows; i++)
        {
            double m = copy.elem[i][k];
            if (i != h)
            {
                for (size_t j = 0; j < copy.cols; j++)
                {
                    copy.elem[i][j] -= m * copy.elem[h][j];
                    inv.elem[i][j] -= m * inv.elem[h][j];
                }
            }
        }
        h++;
        k++;
    }

    matrix id = Id(copy.rows);

    if (matrix_compare(copy, id))
    {
        free_matrix(&id);
        free_matrix(&copy);
        return inv;
    }
    else
    {
        free_matrix(&id);
        free_matrix(&copy);
        fprintf(stderr, "Error: The matrix is not inversible\n");
        exit(EXIT_FAILURE);
    }
}