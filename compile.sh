#gcc -std=c11 main.c -o main.out -lm -L/shome/jonasp/usr/lib/ -lfftw3 -lgsl -lgslcblas
gcc -o main.out main.c -lgsl -lgslcblas -lm -lfftw3 -std=c11 
#memory checks??
#gcc -o main.out main.c -L/home/jonas/Libraries/electric-fence-master/libefence.a -lgsl -lgslcblas -lm -lfftw3 -std=c11 
#-I/usr/include/gsl -L/usr/lib

