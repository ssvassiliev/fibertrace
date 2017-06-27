# icc streamline.c -L/opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64 -llapack -lblas -lgfortran -lpthread -lquadmath -openmp  -o streamline -xHost -static
gcc streamline.c -llapack -lblas -lgfortran -lpthread -lquadmath -lm -fopenmp -O3 -o streamline -Wno-unused-result

