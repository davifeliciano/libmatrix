#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include "matrix.h"

struct eigen
{
    double value;
    double error;
};

struct eigen basic_power(matrix mtx, matrix *x)
{
    matrix prod = matrix_product(mtx, *x);
    double max = max_entry(prod);

    for (size_t i = 0; i < prod.rows; i++)
        prod.elem[i][0] *= 1 / max;

    matrix diff = matrix_diff(prod, *x);
    double error = infinity_norm(diff);

    for (size_t i = 0; i < x->rows; i++)
        x->elem[i][0] = prod.elem[i][0];

    free_matrix(&diff);
    free_matrix(&prod);

    struct eigen result = {max, error};
    return result;
}

double rand_between(double lo, double hi)
{
    return lo + ((double)rand() / (double)RAND_MAX * (hi - lo));
}

double big_eigen(matrix mtx, bool *flag, double (*norm)())
{ /* Return the dominant eigenvalue of a matrix
    If the method converge, set flag as true , false otherwise
    This method is particularly inefficient
    Divergence or slow convergence occurs for a wide class of matrixes */

    if (mtx.rows != mtx.cols)
        error("The first argument in basic_power() must be a square matrix");

    // Creating initial vector
    matrix x = create_matrix(mtx.rows, 1);

    // Initializing initial vector as a random linar combination of cols of mtx
    srand(time(0));
    double n = norm(mtx);

    for (size_t i = 0; i < x.rows; i++)
        x.elem[i][0] = 0.0;

    for (size_t i = 0; i < mtx.rows; i++)
    {
        for (size_t j = 0; j < mtx.cols; j++)
        {
            double c = rand_between(-n / mtx.rows, n / mtx.rows);
            x.elem[i][0] += c * mtx.elem[i][j];
        }
    }

    struct eigen result;
    clock_t tick = clock();
    do
    {
        result = basic_power(mtx, &x);
        // If the computation took more than 4 seconds, exit
        if ((clock() - tick) / CLOCKS_PER_SEC > 4.0)
        {
            *flag = false;
            fprintf(stderr, "Divergence or slow convergence of basic_power method\n");
            return 0.0;
        }
    } while (result.error > THRESHOLD);

    free_matrix(&x);
    *flag = true;
    return result.value;
}

matrix *qr_decompose(matrix mtx)
{
    if (mtx.rows != mtx.cols)
        error("The argument in qr_algorithm() must be a square matrix");

    matrix q = Id(mtx.rows);
    matrix r = copy_matrix(mtx);
    matrix id = Id(mtx.rows);

    for (size_t i = 0; i < r.rows - 1; i++)
    {
        matrix c = get_col(r, i);
        for (size_t j = 0; j < i; j++)
            c.elem[j][0] = 0.0;

        matrix e = create_matrix(r.rows, 1);
        for (size_t j = 0; j < e.rows; j++)
            e.elem[j][0] = 0.0;

        if (c.elem[i][0] > 0)
            e.elem[i][0] = euclidian_norm(c);
        else
            e.elem[i][0] = -euclidian_norm(c);

        matrix v_transp = matrix_sum(c, e);
        matrix v = copy_matrix(v_transp);
        transpose(&v_transp);
        double v_dot_v = dot_product(v, v);
        matrix w = matrix_product(v, v_transp);
        scalar_matrix_product(2 / v_dot_v, &w);

        matrix h = matrix_diff(id, w);
        matrix new_q = matrix_product(q, h);
        matrix new_r = matrix_product(h, r);
        assign_matrix(new_q, &q);
        assign_matrix(new_r, &r);

        free_matrix(&c);
        free_matrix(&e);
        free_matrix(&v);
        free_matrix(&v_transp);
        free_matrix(&w);
        free_matrix(&h);
        free_matrix(&new_q);
        free_matrix(&new_r);
    }

    matrix *result = (matrix *)malloc(2 * sizeof(matrix));
    result[0] = q;
    result[1] = r;

    free_matrix(&id);

    return result;
}

double max_lower_triangle(matrix mtx)
{ // Evaluates the element with greatest absolute value in lower triangle portion of a matrix
    double max = 0;
    for (int j = 0; j < mtx.cols; j++)
    {
        for (int i = j + 1; i < mtx.rows; i++)
        {
            if (fabs(mtx.elem[i][j]) > max)
                max = mtx.elem[i][j];
        }
    }
    return max;
}

double *eigenvalues_qr(matrix mtx)
{
    if (mtx.rows != mtx.cols)
        error("The argument in eigenvalues_qr() must be a square matrix");

    matrix copy = copy_matrix(mtx);

    while (fabs(max_lower_triangle(copy)) > THRESHOLD)
    {
        matrix *qr = qr_decompose(copy);
        matrix prod = matrix_product(qr[1], qr[0]);

        assign_matrix(prod, &copy);
        print_matrix(copy, 12);

        for (size_t i = 0; i < 2; i++)
            free_matrix(&qr[i]);
        free(qr);
        free_matrix(&prod);
    }

    double *result = (double *)malloc(copy.rows * sizeof(double));
    for (size_t i = 0; i < copy.rows; i++)
        result[i] = copy.elem[i][i];

    free_matrix(&copy);

    return result;
}

double *eigenvalues_sm(matrix mtx)
{
    if (mtx.rows != mtx.cols)
        error("The argument of eigenvalues() must be a square matrix");

    matrix id = Id(mtx.rows);
    double *result = (double *)malloc(sizeof(double));

    bool flag1;
    double eigen = big_eigen(mtx, &flag1, infinity_norm);
    if (flag1)
    {
        result[0] = eigen;
        matrix shifted = create_matrix(mtx.rows, mtx.cols);
        size_t count = 1;

        while (count < mtx.rows && flag1)
        {
            for (size_t i = 0; i < id.rows; i++)
                id.elem[i][i] = eigen;

            shifted = matrix_diff(mtx, id);
            eigen = big_eigen(shifted, &flag1, infinity_norm);
            if (flag1)
            {
                result = realloc((void *)result, (count + 1) * sizeof(double));
                result[count] = eigen + result[count - 1];
            }
            else
                printf("Failed on eigenvalue #%ld\n\n", count + 1);
        }

        free_matrix(&shifted);
    }
    else
        fprintf(stderr, "Failed to compute dominant eigenvalue of the given matrix");

    free_matrix(&id);
    return result;
}