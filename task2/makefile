VPATH = ../csv:../:./
CC=gcc

all: task2_eigen task2_qr

matrix.o: matrix.c matrix.h
	$(CC) -c ../matrix.c

libcsv.o: libcsv.c csv.h
	$(CC) -c ../csv/libcsv.c

file_to_matrix.o: file_to_matrix.c matrix.h csv.h
	$(CC) -c ../file_to_matrix.c

eigenvalues.o: matrix.h
	$(CC) -c ../eigenvalues.c

norm.o: matrix.h
	$(CC) -c ../norm.c

task2_qr: matrix.o libcsv.o file_to_matrix.o eigenvalues.o norm.o
	$(CC) matrix.o libcsv.o file_to_matrix.o eigenvalues.o norm.o task2_qr.c -o task2_qr -lm

task2_eigen: matrix.o libcsv.o file_to_matrix.o eigenvalues.o norm.o
	$(CC) matrix.o libcsv.o file_to_matrix.o eigenvalues.o norm.o task2_eigen.c -o task2_eigen -lm

clean:
	rm -rf *.o