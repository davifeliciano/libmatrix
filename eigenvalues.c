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
    double max = 0.0;
    for (int j = 0; j < mtx.cols; j++)
    {
        for (int i = j + 1; i < mtx.rows; i++)
        {
            if (fabs(mtx.elem[i][j]) > max)
                max = fabs(mtx.elem[i][j]);
        }
    }
    return max;
}

matrix eigenvalues_qr(matrix mtx)
{ /* A better algorithm to evaluate eigenvalues
    In the abssense of complex arithmetic, it computes
    only the real eigenvalues
    Double shift QR factorization is the way to go
    if you're trying to compute all eigenvalues using only
    real arithmetic. */

    if (mtx.rows != mtx.cols)
        error("The argument in eigenvalues_qr() must be a square matrix");

    matrix copy = copy_matrix(mtx);
    size_t count = 0;

    while (max_lower_triangle(copy) > THRESHOLD / 2 && count < copy.rows * 1000)
    {
        matrix *qr = qr_decompose(copy);
        matrix prod = matrix_product(qr[1], qr[0]);

        assign_matrix(prod, &copy);

        for (size_t i = 0; i < 2; i++)
            free_matrix(&qr[i]);
        free(qr);
        free_matrix(&prod);

        count++;
    }

    //    printf("shur matrix =\n");
    //    print_matrix(copy, 12);

    size_t result_count = 0; // How many eigenvalues found yet?
    size_t i = 0;
    double *result = (double *)malloc(copy.rows * sizeof(double));

    while (i < copy.rows)
    {
        if (i == copy.rows - 1)
        {
            result[result_count] = copy.elem[i][i];
            result_count++;
            i++;
        }
        else
        {
            if (fabs(copy.elem[i + 1][i]) < THRESHOLD)
            {
                result[result_count] = copy.elem[i][i];
                result_count++;
                i++;
            }
            else
                i += 2;
        }
    }

    if (result_count == 0)
        error("The matrix has no real eigenvalues");

    matrix eigen = create_matrix(result_count, 1);
    for (size_t i = 0; i < eigen.rows; i++)
        eigen.elem[i][0] = result[i];

    free(result);
    free_matrix(&copy);

    return eigen;
}