#include <stddef.h>
#include <stdbool.h>

#ifndef MATRIX_H
#define MATRIX_H

#define THRESHOLD 10e-12 // Any value bellow this one will be seen as a zero

typedef struct
{
    size_t rows;
    size_t cols;
    double **elem;
} matrix;

// Function prototypes
void error(char *msg);
matrix create_matrix(size_t nrow, size_t ncol);
void free_matrix(matrix *mtx);
void print_matrix(matrix mtx, int precision);
matrix copy_matrix(matrix mtx);
void assign_matrix(matrix a, matrix *b);
void transpose(matrix *mtx);
void swap_rows(matrix *mtx, size_t row1, size_t row2);
void swap_cols(matrix *mtx, size_t col1, size_t col2);
int null_row(matrix mtx, size_t row);
int null_col(matrix mtx, size_t col);
void remove_row(matrix *mtx, size_t row);
void remove_col(matrix *mtx, size_t col);
matrix get_row(matrix mtx, size_t row_index);
matrix get_col(matrix mtx, size_t col_index);
matrix Id(size_t n);
int matrix_compare(matrix a, matrix b);
matrix matrix_sum(matrix a, matrix b);
matrix matrix_diff(matrix a, matrix b);
matrix matrix_product(matrix a, matrix b);
void scalar_matrix_product(double x, matrix *mtx);
double dot_product(matrix a, matrix b);

void gauss_eliminate(matrix *mtx, size_t args, ...);
matrix solve(matrix *mtx, matrix *b);
double det(matrix mtx);

matrix inverse(matrix mtx);

matrix import_data(char filename[]);
void matrix_to_file(matrix mtx, int precision, char *filename);

matrix *lu_decompose(matrix mtx);
matrix solve_lower(matrix a, matrix b);
matrix solve_upper(matrix a, matrix b);

double max_entry(matrix mtx);
double infinity_norm(matrix mtx);
double euclidian_norm(matrix mtx);

double big_eigen(matrix mtx, bool *flag, double (*norm)());
double small_eigen(matrix mtx, bool *flag, double (*norm)());
matrix *qr_decompose(matrix mtx);
double *eigenvalues_qr(matrix mtx);

#endif