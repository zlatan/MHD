gcc -Wall -I/usr/local/include -c matrix.c 
gcc -L/usr/local/lib matrix.o -lgsl -lgslcblas -lm
