#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "../matrix.h"

float random_float(float min, float max)
{
    return (float)rand() / RAND_MAX * (max - min) + min;
}

void print_usage(char *program_name)
{
    fprintf(stderr, "Usage: %s [-s SEED] [-o OUTPUT] [ROWS] [COLS]\n", program_name);
    fprintf(stderr, "Options are optional\n");
    exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{ /* Generate a random matrix
    Dimensions are passed as 2 integer cl arguments
    You can provide some integer seed with -s option
    You can provide the output file with -o option */
    char ch;
    int seed = 0;
    char *filename = NULL;
    char *program = (char *)malloc((strlen(argv[0]) + 1) * sizeof(char));
    strcpy(program, argv[0]);

    while ((ch = getopt(argc, argv, "s:o:")) != EOF)
    {
        switch (ch)
        {
        case 's':
            seed = abs(atoi(optarg));
            break;

        case 'o':
            filename = optarg;
            break;

        default:
            print_usage(program);
            break;
        }
    }

    argc -= optind;
    argv += optind;

    if (argc < 2)
    {
        fprintf(stderr, "Too few arguments\n\n");
        print_usage(program);
    }

    if (argc > 2)
    {
        fprintf(stderr, "Too many arguments\n\n");
        print_usage(program);
    }

    if (seed)
        srand(seed);
    else
        srand(time(0));

    int dims[2] = {0, 0};
    for (size_t i = 0; i < 2; i++)
    {
        dims[i] = atoi(argv[i]);
        if (dims[i] <= 0)
        {
            fprintf(stderr, "The arguments must be positive integers\n\n");
            print_usage(program);
        }
    }

    matrix a = create_matrix(dims[0], dims[1]);

    for (size_t i = 0; i < a.rows; i++)
    {
        for (size_t j = 0; j < a.cols; j++)
            a.elem[i][j] = random_float(-10, 10);
    }

    if (filename)
    {
        matrix_to_file(a, 5, filename);
        if (seed)
            printf("filename = %s; seed = %d\n", filename, seed);
        else
            printf("filename = %s\n", filename);
    }
    else
    {
        if (seed)
            printf("seed = %d\n", seed);
    }

    print_matrix(a, 5);
    return 0;
}