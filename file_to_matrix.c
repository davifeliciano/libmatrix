#include <stdio.h>
#include <stdlib.h>
#include "csv/csv.h"
#include "matrix.h"

struct counts
{
    char *filename;
    size_t rows;
    size_t cols_count;
    size_t cols;
};

struct foo
{
    matrix mtx;
    struct counts count;
};

void first_cb1(void *s, size_t len, void *data)
{
    ((struct counts *)data)->cols_count++;
}

void first_cb2(int c, void *data)
{
    ((struct counts *)data)->rows++;

    if (((struct counts *)data)->cols)
    {

        if (((struct counts *)data)->cols != ((struct counts *)data)->cols_count)
        {
            fprintf(stderr, "Error: The file %s is not a matrix\n", ((struct counts *)data)->filename);
            exit(EXIT_FAILURE);
        }
        else
        {
            ((struct counts *)data)->cols = ((struct counts *)data)->cols_count;
            ((struct counts *)data)->cols_count = 0;
        }
    }
    else
    {
        ((struct counts *)data)->cols = ((struct counts *)data)->cols_count;
        ((struct counts *)data)->cols_count = 0;
    }
}

void second_cb1(void *s, size_t len, void *data)
{

    size_t current_row = ((struct foo *)data)->count.rows;
    size_t current_col = ((struct foo *)data)->count.cols_count;
    ((struct foo *)data)->mtx.elem[current_row][current_col] = atof((char *)s);
    ((struct foo *)data)->count.cols_count++;
}

void second_cb2(int c, void *data)
{
    ((struct foo *)data)->count.rows++;

    if (((struct foo *)data)->count.cols)
    {
        ((struct foo *)data)->count.cols = ((struct foo *)data)->count.cols_count;
        ((struct foo *)data)->count.cols_count = 0;
    }
    else
    {
        ((struct foo *)data)->count.cols = ((struct foo *)data)->count.cols_count;
        ((struct foo *)data)->count.cols_count = 0;
    }
}

matrix import_data(char filename[])
{
    FILE *input;
    struct csv_parser p;
    char buffer[1024];
    size_t bytes_read;
    struct counts c = {filename, 0, 0, 0};

    if (csv_init(&p, CSV_APPEND_NULL))
    {
        fprintf(stderr, "Error: Failed to initialize csv parser\n");
        exit(EXIT_FAILURE);
    }

    input = fopen(filename, "rb");

    if (!input)
    {
        fprintf(stderr, "Error: Failed to open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    while ((bytes_read = fread(buffer, 1, 1024, input)) > 0)
    {
        if (csv_parse(&p, buffer, bytes_read, first_cb1, first_cb2, &c) != bytes_read)
        {
            fprintf(stderr, "Error: Failed to parse the file %s: %s\n",
                    filename, csv_strerror(csv_error(&p)));
        }
    }

    csv_fini(&p, first_cb1, first_cb2, &c);
    csv_free(&p);
    fclose(input);

    matrix m = create_matrix(c.rows, c.cols);

    c.rows = 0;
    c.cols = 0;
    c.cols_count = 0;

    struct foo foobar = {m, c};

    if (csv_init(&p, 0))
    {
        fprintf(stderr, "Error: Failed to initialize csv parser\n");
        exit(EXIT_FAILURE);
    }

    input = fopen(filename, "rb");

    if (!input)
    {
        fprintf(stderr, "Error: Failed to open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    while ((bytes_read = fread(buffer, 1, 1024, input)) > 0)
    {
        if (csv_parse(&p, buffer, bytes_read, second_cb1, second_cb2, &foobar) != bytes_read)
        {
            fprintf(stderr, "Error: Failed to parse the file %s: %s\n",
                    filename, csv_strerror(csv_error(&p)));
        }
    }

    csv_fini(&p, second_cb1, second_cb2, &foobar);
    csv_free(&p);

    fclose(input);

    return foobar.mtx;
}