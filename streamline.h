#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/timeb.h>

#pragma warning(disable:981)

float*** gridXX;
float*** gridXY;
float*** gridYY;
float*** gridXZ; 
float*** gridYZ; 
float*** gridZZ;

int nb_elementsX, nb_elementsY, nb_elementsZ;
float x_min, y_min, z_min;
float res, dt; 

// FILTERING
float sigma; // Half width of the filter at 1/e
int nb_points; // Number of intergration intervals 

typedef struct {
        double c[3][3];
	int flag;
} TENSOR;

TENSOR*** fields;

typedef struct point {
  float x; float y; float z;
  float a; float c; float b;
} P;

typedef struct {
        float egv_x[3], egv_y[3], egv_z[3];
        float e[3]; float ADC; float FA;
} EIGEN;

typedef struct {
        int *val;
        int size;
} matrix;

P* trace, *trace_tmp;

//Eigenvalues & orthonormal eigenvectors of symmetric matrix fortran code
extern void dsyev_();
//for linear solving
extern void sgetrf_(); extern void sgetrs_(); extern void sppsv_();extern void dppsv_();extern void dposv_();
// Map xyz directions to float using 1000 color indexes 
double xyzmap1000(float *);
 // Compute eigenvalues & eigenvectors of diffusion tensor using LAPACK lib
EIGEN eig_lapack(TENSOR);
// Make sample tensor fiels for debugging
void init_sample(void);
// Read tensor field
void read_tensor_file(char*);
// Compute tensor at point po  using trilinear interpolation
TENSOR interpl_tens(P);
TENSOR mlsfilter(EIGEN, P, int, float, float, int);
inline float mountain(EIGEN);
inline float* direction(EIGEN); 
inline float dot_product(float *, float *);
inline void invert_vector(float *);
void save_interp_grid(int);
void save_streamline_pdb(P*, int *, long int);
void save_box(int, int, int, int, int, int);
