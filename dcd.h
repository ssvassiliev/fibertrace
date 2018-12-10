typedef struct H_ {
        int charmm,charmm_extrablock,charmm_4dims,NSet,ISTART,NSAVC,NAMNF,N,FREEINDEXES,endoffile;
        FILE* fp;
        float DELTA;
} H;

int read_dcdheader(char* filename);
int write_crdfile(char* filename);
int read_dcdstep(float *x, float *y, float *z);

int read_xyz(int frame_num, float *x, float *y, float*z);
int reset_file();
int initialize_read_xyz(char* filename);
int forward_dcdstep();

int close_dcd();
