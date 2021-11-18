VPATH = ../csv:../:./

all: lu ge

matrix.o: matrix.c matrix.h
	gcc -c ../matrix.c

libcsv.o: libcsv.c csv.h
	gcc -c ../csv/libcsv.c

file_to_matrix.o: file_to_matrix.c matrix.h csv.h
	gcc -c ../file_to_matrix.c

gauss_elimination.o: gauss_elimination.c matrix.h
	gcc -c ../gauss_elimination.c

lu_decomposition.o: lu_decomposition.c matrix.h
	gcc -c ../lu_decomposition.c

lu: matrix.o libcsv.o file_to_matrix.o lu_decomposition.o
	gcc matrix.o libcsv.o file_to_matrix.o lu_decomposition.o lu.c -o lu

ge: matrix.o libcsv.o file_to_matrix.o gauss_elimination.o
	gcc matrix.o libcsv.o file_to_matrix.o gauss_elimination.o ge.c -o ge

clean:
	rm -rf *.o