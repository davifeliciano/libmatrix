#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

double max_entry(matrix mtx)
{
    double max = mtx.elem[0][0];
    for (size_t i = 1; i < mtx.rows; i++)
    {
        for (size_t j = 0; j < mtx.cols; j++)
        {
            if (mtx.elem[i][0] > max)
                max = mtx.elem[i][j];
        }
    }

    return max;
}

double infinity_norm(matrix mtx)
{
    matrix rows_sum = create_matrix(mtx.rows, 1);
    double norm;

    for (size_t i = 0; i < mtx.rows; i++)
    {
        double sum = 0;
        for (size_t j = 0; j < mtx.cols; j++)
            sum += fabs(mtx.elem[i][j]);
        rows_sum.elem[i][0] = sum;
    }

    norm = max_entry(rows_sum);
    free_matrix(&rows_sum);
    return norm;
}