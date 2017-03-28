typedef struct {
        float bx,by,bz;
        float xmin,xmax;
        float ymin,ymax;
        float zmin,zmax;
        float density;
        float timescale;
        int fstep;
        float cutoff;
        char** dcd_files;
        int files_count;
        char* prmtop;
        char adc;
        char diff;
        char fa;
        char hydro;
        char oxy;
        char tensors;
} Setup;

typedef struct {
        float xx; //float yx; float zx;
        float xy; float yy; //float zy;
        float xz; float yz; float zz;
} DISP;

typedef struct {
	float egv_x[3], egv_y[3], egv_z[3];
	float e[3];
	float ADC;
	float FA;
	int b;
} EIGEN;

inline float Dab(float v1_t1, float v1_t2, float v2_t1, float v2_t2, float T);
void WaterDensity(Setup* S);
inline void transpose_square(double** I);
