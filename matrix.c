#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

matrix create_matrix(size_t nrow, size_t ncol) {
	matrix mtx;
    mtx.rows = nrow;
    mtx.cols = ncol;
	mtx.elem = (double **) malloc((size_t) (nrow * sizeof(double*)));
    mtx.elem[0] = (double *) malloc((size_t) (nrow * ncol * sizeof(double)));    
    for (size_t i = 1; i < nrow; i++)
        mtx.elem[i] = mtx.elem[i - 1] + ncol;
    return mtx;
}

void free_matrix(matrix *mtx) {
    free((void *) (mtx->elem[0]));
    free((void *) mtx->elem);
} 

void print_matrix(matrix mtx, int precision) {
    for (size_t i = 0; i < mtx.rows; i++) {
        for (size_t j = 0; j < mtx.cols; j++) {
            if (j != mtx.cols - 1)
                printf("%.*f\t", precision, mtx.elem[i][j]);
            else
                printf("%.*f\n", precision, mtx.elem[i][j]);
        }
    } 
}

matrix copy_matrix(matrix mtx) {
    matrix copy = create_matrix(mtx.rows, mtx.cols);
    for (size_t i = 0; i < copy.rows; i++) {
        for (size_t j = 0; j < copy.cols; j++)
            copy.elem[i][j] = mtx.elem[i][j];
    }
    return copy;    
}

void transpose(matrix *mtx) {
    matrix tmp_mtx = copy_matrix(*mtx);
    mtx->rows = tmp_mtx.cols;
    mtx->cols = tmp_mtx.rows;
    mtx->elem = realloc((void *)(mtx->elem), (size_t) (mtx->rows * sizeof(double*)));
    for (size_t i = 1; i < mtx->rows; i++)
        mtx->elem[i] = mtx->elem[i - 1] + mtx->cols;
    for (size_t i = 0; i < mtx->rows; i++) {
        for (size_t j = 0; j < mtx->cols; j++)
            mtx->elem[i][j] = tmp_mtx.elem[j][i];
    }
    free_matrix(&tmp_mtx);
}

void swap_rows(matrix *mtx, size_t row1, size_t row2) {
    if (row1 >= mtx->rows || row2 >= mtx->rows) {
        fprintf(stderr, "Error: Argument exceeds the dimension of the matrix in swap_rows().\n");
        exit(EXIT_FAILURE);
    }
    double tmp_row[mtx->cols];
    for (size_t i = 0; i < mtx->cols; i++) {
        tmp_row[i] = mtx->elem[row1][i];
        mtx->elem[row1][i] = mtx->elem[row2][i];
        mtx->elem[row2][i] = tmp_row[i];
    }
}

void swap_cols(matrix *mtx, size_t col1, size_t col2) {
    if (col1 >= mtx->cols || col2 >= mtx->cols) {
        fprintf(stderr, "Error: Argument exceeds the dimension of the matrix in swap_cols()\n");
        exit(EXIT_FAILURE);
    }
    double tmp_col[mtx->rows];
    for (size_t i = 0; i < mtx->rows; i++) {
        tmp_col[i] = mtx->elem[i][col1];
        mtx->elem[i][col1] = mtx->elem[i][col2];
        mtx->elem[i][col2] = tmp_col[i];
    }
}  

int null_row(matrix mtx, size_t row) {
    for (size_t j = 0; j < mtx.cols ; j++) {
        if (mtx.elem[row][j])
            return 0;
    }
    return 1;
}

int null_col(matrix mtx, size_t col) {
    for (size_t i = 0; i < mtx.rows ; i++) {
        if (mtx.elem[i][col])
            return 0;
    }
    return 1;
}

matrix get_row(matrix mtx, size_t row_index) {
    if (row_index >= mtx.rows) {
        fprintf(stderr, "Error: Argument exceeds the dimension of the matrix in get_row().\n");
        exit(EXIT_FAILURE);
    }
    matrix row = create_matrix(1, mtx.cols);
    for (size_t i = 0; i < mtx.cols; i++)
        row.elem[0][i] = mtx.elem[row_index][i];
    return row;
}

matrix get_col(matrix mtx, size_t col_index) {
    if (col_index >= mtx.cols) {
        fprintf(stderr, "Error: Argument exceeds the dimension of the matrix in get_col()\n");
        exit(EXIT_FAILURE);       
    }
    matrix col = create_matrix(mtx.rows, 1);
    for (size_t i = 0; i < mtx.rows; i++)
        col.elem[i][0] = mtx.elem[i][col_index];
    return col;
}

void remove_row(matrix *mtx, size_t row) {
    if (row > mtx->rows) {
        fprintf(stderr, "Error: Argument exceeds the dimension of the matrix in remove_row()\n");
        exit(EXIT_FAILURE);
    }
    if (mtx->rows == 1) {
        fprintf(stderr, "Error: A matriz tem apenas uma linha in remove_row()\n");
        exit(EXIT_FAILURE);
    }
    matrix copy = copy_matrix(*mtx);
    mtx->rows--;
    mtx->elem = realloc((void *)(mtx->elem), (size_t) (mtx->rows * sizeof(double*)));
    for (size_t i = 1; i < mtx->rows; i++)
        mtx->elem[i] = mtx->elem[i - 1] + mtx->cols;
    for (size_t i = 0; i < mtx->rows; i++) {
        for (size_t j = 0; j < mtx->cols; j++) {
            if (i < row)
                mtx->elem[i][j] = copy.elem[i][j];    
            else
                mtx->elem[i][j] = copy.elem[i + 1][j];
        }
    }
    free_matrix(&copy);
}

void remove_col(matrix *mtx, size_t col) {
    if (col > mtx->cols) {
        fprintf(stderr, "Error: Argument exceeds the dimension of the matrix in remove_col()\n");
        exit(EXIT_FAILURE);
    }
    if (mtx->cols == 1) {
        fprintf(stderr, "Error: The matrix has just one column in remove_col()\n");
        exit(EXIT_FAILURE);
    }
    matrix copy = copy_matrix(*mtx);
    mtx->cols--;
    mtx->elem[0] = realloc((void *)(mtx->elem[0]), (size_t) (mtx->rows * mtx->cols * sizeof(double)));
    for (size_t i = 1; i < mtx->rows; i++)
        mtx->elem[i] = mtx->elem[i - 1] + mtx->cols;
    for (size_t i = 0; i < mtx->rows; i++) {
        for (size_t j = 0; j < mtx->cols; j++) {
            if (j < col)
                mtx->elem[i][j] = copy.elem[i][j];    
            else
                mtx->elem[i][j] = copy.elem[i][j + 1];
        }
    }
    free_matrix(&copy);
}

matrix Id(size_t n) {
    matrix id = create_matrix(n, n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (i == j)
                id.elem[i][j] = 1;
            else
                id.elem[i][j] = 0;
        }   
    }
    return id;
}

int matrix_compare(matrix a, matrix b) {
    if (a.rows != b.rows)
        return 0;
    if (a.cols != b.cols)
        return 0;
    for (size_t i = 0; i < a.rows; i++) {
        for (size_t j = 0; j < b.cols; j++) {
            if (fabs(a.elem[i][j] - b.elem[i][j]) > 1.0E-12)
                return 0;
        }
    }
    return 1;
}

matrix matrix_sum(matrix a, matrix b) {
    if (a.rows != b.rows || a.cols != b.cols) {
        fprintf(stderr, "Error: The matrixes has different dimensions in matrix_sum()\n");
        exit(EXIT_FAILURE);
    }
    matrix sum = create_matrix(a.rows, a.cols);
    for (size_t i = 0; i < a.rows; i++) {
        for (size_t j = 0; j < a.cols; j++) {
            sum.elem[i][j] = a.elem[i][j] + b.elem[i][j];
        }
    }
    return sum;
}

matrix matrix_diff(matrix a, matrix b) {
    if (a.rows != b.rows || a.cols != b.cols) {
        fprintf(stderr, "Error: The matrixes has different dimensions in matrix_diff()\n");
        exit(EXIT_FAILURE);
    }
    matrix diff = create_matrix(a.rows, a.cols);
    for (size_t i = 0; i < a.rows; i++) {
        for (size_t j = 0; j < a.cols; j++) {
            diff.elem[i][j] = a.elem[i][j] - b.elem[i][j];
        }
    }
    return diff;
}

matrix matrix_product(matrix a, matrix b) {
    if (a.cols != b.rows) {
        fprintf(stderr, "Error: The product is not well defined given the dimensions of the matrix in matrix_product()\n");
        exit(EXIT_FAILURE);
    }
    matrix prod = create_matrix(a.rows, b.cols);
    for (size_t i = 0; i < prod.rows; i++) {
        for (size_t j = 0; j < prod.cols; j++) {
            float sum = 0;
            for (size_t k = 0; k < a.cols; k++)
                sum += a.elem[i][k] * b.elem[k][j];
            prod.elem[i][j] = sum;
        }
    }
    return prod;
}