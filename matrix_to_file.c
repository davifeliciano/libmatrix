#include <stdio.h>
#include <stdlib.h>
#include "csv/csv.h"
#include "matrix.h"

int count_digits(double n)
{
    if (abs(n) < 1)
        return 1;
    int count = 0;
    while (abs(n) >= 1)
    {
        n = n / 10;
        count++;
    }
    return count;
}

void matrix_to_file(matrix mtx, int precision, char *filename)
{ /* Export a given matrix to a file in csv format, provided the filename
    The number of decimal cases is defined by precision argument */
    FILE *output;
    output = fopen(filename, "wb");
    if (!output)
    {
        fprintf(stderr, "Error: Failed to open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    char buffer[64];
    double elem;
    int digits;

    for (size_t i = 0; i < mtx.rows; i++)
    {
        for (size_t j = 0; j < mtx.cols - 1; j++)
        {
            elem = mtx.elem[i][j];
            digits = count_digits(elem);
            sprintf(buffer, "%.*f", precision, elem);

            if (csv_fwrite(output, (void *)buffer, precision + 1 + ((elem < 0) ? digits + 1 : digits)) == EOF)
            // decimal digits + radix char + digits + minus if negative
            {
                fprintf(stderr, "Error while writing %s", filename);
                exit(EXIT_FAILURE);
            }
            fputc(44, output); // put a coma
        }

        elem = mtx.elem[i][mtx.cols - 1];
        digits = count_digits(elem);
        sprintf(buffer, "%.*f", precision, elem);

        if (csv_fwrite(output, (void *)buffer, precision + 1 + ((elem < 0) ? digits + 1 : digits)) == EOF)
        {
            fprintf(stderr, "Error while writing %s", filename);
            exit(EXIT_FAILURE);
        }
        fputc(13, output); // put a return char
    }
    fclose(output);
}
