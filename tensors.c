#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>

#include "dcd.c"
#include "file_read_write.h"
#include "file_read_write.c"
#include "tensors.h"

extern void dsyev_();//eigen values & orthonormal eigenvectors* fortran code

inline static float Dab(float v1_t1, float v1_t2, float v2_t1, float v2_t2, float T) {
  return (( (v1_t2-v1_t1) * (v2_t2-v2_t1) ) / (2.0f*T));}

void WaterDensity(Setup* S) {

  int i,j,o,k,m;
  //parameters for LAPACK
  char JOBZ = 'V';
  char UPLO = 'L';
  int N = 3;
  double AA[3][3];
  int LDA = 3;
  double W[3];
  double WORK[4*3];
  int LWORK = 4*3;
  int INFO;
  float res = S->density;
  float timescale = S->timescale;
  float cut_off = S->cutoff;

  if(read_parm7(S->prmtop)==1)
    exit(1);
  if(!initialize_read_xyz(S->dcd_files[0])) return; // dcd_file
  FILE* fp1;
  FILE* fp2;
  FILE* fp3;
  FILE* fp4;
  FILE* fp5;
  int frames = my_H.NSet;
  int nb_atoms = my_H.N;
  int nb_elementsX,nb_elementsY,nb_elementsZ;
  int wat,x_index,y_index,z_index;
  float X[2][nb_atoms]; float Y[2][nb_atoms]; float Z[2][nb_atoms];

  //cube dimensions
  float x_min = S->xmin; float x_max = S->xmax;
  float y_min = S->ymin; float y_max = S->ymax;
  float z_min = S->zmin; float z_max = S->zmax;

  //gathering the water indices
  int water_count = 0;
  int* water_ptr = calloc(0,sizeof(int));
  for(j=0; j<nb_atoms; j++) {
    if (!strncmp("WAT",resName[j],3) && !strncmp("O",atom_name[j],1)) {
      water_ptr = realloc(water_ptr,sizeof(int)*(water_count+1));
      water_ptr[water_count] = j; water_count++;
    }
  }
  printf("Number of Oxygens: %i\n",water_count);
  printf("Number of Hydrogens: %i\n",water_count*2);

  nb_elementsX = (int)((x_max-x_min)/S->density); nb_elementsY = (int)((y_max-y_min)/S->density); nb_elementsZ = (int)((z_max-z_min)/S->density);
  int***  tensor_count;

  float*** gridDS;
  float*** grid;
  float*** gridO;
  float*** gridH;

  double*** gridXX;
  double*** gridXY;
  double*** gridYY;
  double*** gridXZ;
  double*** gridYZ;
  double*** gridZZ;
  tensor_count = (int***)calloc(nb_elementsX ,sizeof(int**));
  gridDS = (float***)calloc(nb_elementsX ,sizeof(float**));
  grid   = (float***)calloc(nb_elementsX ,sizeof(float**));
  gridO  = (float***)calloc(nb_elementsX ,sizeof(float**));
  gridH  = (float***)calloc(nb_elementsX ,sizeof(float**));

  gridXX = (double***)calloc(nb_elementsX ,sizeof(double**));
  gridXY = (double***)calloc(nb_elementsX ,sizeof(double**));
  gridYY = (double***)calloc(nb_elementsX ,sizeof(double**));
  gridXZ = (double***)calloc(nb_elementsX ,sizeof(double**));
  gridYZ = (double***)calloc(nb_elementsX ,sizeof(double**));
  gridZZ = (double***)calloc(nb_elementsX ,sizeof(double**));
  if (tensor_count==NULL) {printf("too big!\n"); exit(0);}
  if (gridDS==NULL) {printf("too big!\n"); exit(0);}
  if (grid==NULL)   {printf("too big!\n"); exit(0);}
  if (gridO==NULL)  {printf("too big!\n"); exit(0);}
  if (gridH==NULL)  {printf("too big!\n"); exit(0);}
  if (gridXX==NULL) {printf("too big!\n"); exit(0);}
  if (gridXY==NULL) {printf("too big!\n"); exit(0);}
  if (gridYY==NULL) {printf("too big!\n"); exit(0);}
  if (gridXZ==NULL) {printf("too big!\n"); exit(0);}
  if (gridYZ==NULL) {printf("too big!\n"); exit(0);}
  if (gridZZ==NULL) {printf("too big!\n"); exit(0);}
  for(i=0; i<nb_elementsX; i++){
    tensor_count[i] = (int**)calloc(nb_elementsY,sizeof(int*));
    gridDS[i] = (float**)calloc(nb_elementsY,sizeof(float*));
    grid[i]   = (float**)calloc(nb_elementsY,sizeof(float*));
    gridO[i]  = (float**)calloc(nb_elementsY,sizeof(float*));
    gridH[i]  = (float**)calloc(nb_elementsY,sizeof(float*));
    gridXX[i] = (double**)calloc(nb_elementsY,sizeof(double*));
    gridXY[i] = (double**)calloc(nb_elementsY,sizeof(double*));
    gridYY[i] = (double**)calloc(nb_elementsY,sizeof(double*));
    gridXZ[i] = (double**)calloc(nb_elementsY,sizeof(double*));
    gridYZ[i] = (double**)calloc(nb_elementsY,sizeof(double*));
    gridZZ[i] = (double**)calloc(nb_elementsY,sizeof(double*));
    if (tensor_count[i]==NULL) {printf("too big!\n"); exit(0);}
    if (gridDS[i]==NULL) {printf("too big!\n"); exit(0);}
    if (grid[i]==NULL)   {printf("too big!\n"); exit(0);}
    if (gridO[i]==NULL)  {printf("too big!\n"); exit(0);}
    if (gridH[i]==NULL)  {printf("too big!\n"); exit(0);}
    if (gridXX[i]==NULL) {printf("too big!\n"); exit(0);}
    if (gridXY[i]==NULL) {printf("too big!\n"); exit(0);}
    if (gridYY[i]==NULL) {printf("too big!\n"); exit(0);}
    if (gridXZ[i]==NULL) {printf("too big!\n"); exit(0);}
    if (gridYZ[i]==NULL) {printf("too big!\n"); exit(0);}
    if (gridZZ[i]==NULL) {printf("too big!\n"); exit(0);}
    for(j=0; j<nb_elementsY; j++) {
      tensor_count[i][j] = (int*)calloc(nb_elementsZ,sizeof(int));
      gridDS[i][j] = (float*)calloc(nb_elementsZ,sizeof(float));
      grid[i][j]   = (float*)calloc(nb_elementsZ,sizeof(float));
      gridO[i][j]  = (float*)calloc(nb_elementsZ,sizeof(float));
      gridH[i][j]  = (float*)calloc(nb_elementsZ,sizeof(float));
      gridXX[i][j] = (double*)calloc(nb_elementsZ,sizeof(double));
      gridXY[i][j] = (double*)calloc(nb_elementsZ,sizeof(double));
      gridYY[i][j] = (double*)calloc(nb_elementsZ,sizeof(double));
      gridXZ[i][j] = (double*)calloc(nb_elementsZ,sizeof(double));
      gridYZ[i][j] = (double*)calloc(nb_elementsZ,sizeof(double));
      gridZZ[i][j] = (double*)calloc(nb_elementsZ,sizeof(double));
      if (tensor_count[i][j]==NULL) {printf("too big!\n"); exit(0);}
      if (gridDS[i][j]==NULL) {printf("too big!\n"); exit(0);}
      if (grid[i][j]==NULL)   {printf("too big!\n"); exit(0);}
      if (gridO[i][j]==NULL)  {printf("too big!\n"); exit(0);}
      if (gridH[i][j]==NULL)  {printf("too big!\n"); exit(0);}
      if (gridXX[i][j]==NULL) {printf("too big!\n"); exit(0);}
      if (gridXY[i][j]==NULL) {printf("too big!\n"); exit(0);}
      if (gridYY[i][j]==NULL) {printf("too big!\n"); exit(0);}
      if (gridXZ[i][j]==NULL) {printf("too big!\n"); exit(0);}
      if (gridYZ[i][j]==NULL) {printf("too big!\n"); exit(0);}
      if (gridZZ[i][j]==NULL) {printf("too big!\n"); exit(0);}
    }
  }

  for(i=0; i<nb_elementsX; i++)
    for(j=0; j<nb_elementsY; j++)
      for(k=0; k<nb_elementsZ; k++){
	grid[i][j][k] = 0.0;
	gridDS[i][j][k] = gridO[i][j][k] = gridH[i][j][k]=0.0;
	gridXX[i][j][k] = gridXY[i][j][k] = gridYY[i][j][k]=0.0;
	gridXZ[i][j][k] = gridYZ[i][j][k] = gridZZ[i][j][k]=0.0f;
      }

  float total_frames = 0.0f;
  DISP displacement;
  int t1_index = 0; int t2_index = 0;
  float x_avg[2], y_avg[2], z_avg[2];
  float maxDS = -1.0f;
  int idx=0;

  printf("\n\nComputing water density..\n");
  for(m=0; m<S->files_count; m++)
    {
      printf("File %i/%i, Frames: %i\n", m+1, S->files_count, frames);
      for(i=0; i<frames; i++)
	{
	  read_xyz(i,X[0],Y[0],Z[0]);
	  printf("Processing frame %i\r",i);
	  for(j=0; j<water_count; j++)
	    {
	      wat = water_ptr[j];
	      //computing for diffusion output
	      x_avg[0] = (X[0][wat] + X[0][wat+1] + X[0][wat+2]) / 3.0f;
	      y_avg[0] = (Y[0][wat] + Y[0][wat+1] + Y[0][wat+2]) / 3.0f;
	      z_avg[0] = (Z[0][wat] + Z[0][wat+1] + Z[0][wat+2]) / 3.0f;
	      if(x_avg[0]>=x_min && x_avg[0]<=x_max)
		if(y_avg[0]>=y_min && y_avg[0]<=y_max)
		  if(z_avg[0]>=z_min && z_avg[0]<=z_max) {
		    x_index = (int) nearbyint((double)((x_avg[0]-x_min)/res)); if (x_index>=nb_elementsX) x_index = nb_elementsX-1;
		    y_index = (int) nearbyint((double)((y_avg[0]-y_min)/res)); if (y_index>=nb_elementsY) y_index = nb_elementsY-1;
		    z_index = (int) nearbyint((double)((z_avg[0]-z_min)/res)); if (z_index>=nb_elementsZ) z_index = nb_elementsZ-1;
		    gridDS[x_index][y_index][z_index] += 1.0f;
		    if (gridDS[x_index][y_index][z_index] > maxDS) maxDS = gridDS[x_index][y_index][z_index];
		  }
	      //computing for density output
	      //getting oxygen
	      if(X[0][wat]>=x_min && X[0][wat]<=x_max)
		if(Y[0][wat]>=y_min && Y[0][wat]<=y_max)
		  if(Z[0][wat]>=z_min && Z[0][wat]<=z_max) {
		    x_index = nearbyintf((X[0][wat]-x_min)/res); if (x_index>=nb_elementsX) x_index = nb_elementsX-1;
		    y_index = nearbyintf((Y[0][wat]-y_min)/res); if (y_index>=nb_elementsY) y_index = nb_elementsY-1;
		    z_index = nearbyintf((Z[0][wat]-z_min)/res); if (z_index>=nb_elementsZ) z_index = nb_elementsZ-1;
		    gridO[x_index][y_index][z_index] +=1.0f;
		  }
	      //getting hydrogen1
	      if(X[0][wat+1]>=x_min && X[0][wat+1]<=x_max)
		if(Y[0][wat+1]>=y_min && Y[0][wat+1]<=y_max)
		  if(Z[0][wat+1]>=z_min && Z[0][wat+1]<=z_max) {
		    x_index = nearbyintf((X[0][wat+1]-x_min)/res); if (x_index>=nb_elementsX) x_index = nb_elementsX-1;
		    y_index = nearbyintf((Y[0][wat+1]-y_min)/res); if (y_index>=nb_elementsY) y_index = nb_elementsY-1;
		    z_index = nearbyintf((Z[0][wat+1]-z_min)/res); if (z_index>=nb_elementsZ) z_index = nb_elementsZ-1;
		    gridH[x_index][y_index][z_index] +=1.0f;
		  }
	      //getting hydrogen2
	      if(X[0][wat+2]>=x_min && X[0][wat+2]<=x_max)
		if(Y[0][wat+2]>=y_min && Y[0][wat+2]<=y_max)
		  if(Z[0][wat+2]>=z_min && Z[0][wat+2]<=z_max) {
		    x_index = nearbyintf((X[0][wat+2]-x_min)/res); if (x_index>=nb_elementsX) x_index = nb_elementsX-1;
		    y_index = nearbyintf((Y[0][wat+2]-y_min)/res); if (y_index>=nb_elementsY) y_index = nb_elementsY-1;
		    z_index = nearbyintf((Z[0][wat+2]-z_min)/res); if (z_index>=nb_elementsZ) z_index = nb_elementsZ-1;
		    gridH[x_index][y_index][z_index] +=1.0f;
		  }
	    } // waters
	} // frames
      printf("\n");
      close_dcd();
      idx = m+1; if (idx<S->files_count) {
      	initialize_read_xyz(S->dcd_files[idx]);
      	frames = my_H.NSet;
      }
    } // files


  for(i=0; i<nb_elementsX; i++)
    for(j=0; j<nb_elementsY; j++)
      for(k=0; k<nb_elementsZ; k++) {
       	gridDS[i][j][k] /= maxDS;
	tensor_count[i][j][k] = 0.0;
      }
  printf("\n\nDoing diffusion..\n");


  // --------- computing tensor field --------------------
  int scale = S->fstep;//minimum =1;
  int fc;
  int globalFrameCounter=0;
  for(m=0; m< S->files_count; m++) {
    initialize_read_xyz(S->dcd_files[m]);
    frames = my_H.NSet;
    total_frames += (frames/scale);

    printf("File %i/%i, Frames: %i\n", m+1, S->files_count, frames);
    int ij=0;
    for(fc=0;fc<scale;fc++)
      {
	read_xyz(fc,X[t1_index],Y[t1_index],Z[t1_index]);
	for(i=scale+fc; i<frames-1-fc; i+=scale) {
	  t2_index = t1_index^1;
	  read_xyz(i,X[t2_index],Y[t2_index],Z[t2_index]);
	  printf("Processing frame %i, GlobalFrameCounter=%i\r",ij++,globalFrameCounter);
          globalFrameCounter++;
	  for(j=0; j<water_count; j++) {
	    wat = water_ptr[j];
	    // Find center water molecule - size weghted
	    x_avg[t1_index] = (X[t1_index][wat] + 0.5*(X[t1_index][wat+1] + X[t1_index][wat+2]))*0.5;
	    y_avg[t1_index] = (Y[t1_index][wat] + 0.5*(Y[t1_index][wat+1] + Y[t1_index][wat+2]))*0.5 ;
	    z_avg[t1_index] = (Z[t1_index][wat] + 0.5*(Z[t1_index][wat+1] + Z[t1_index][wat+2]))*0.5;
	    x_avg[t2_index] = (X[t2_index][wat] + 0.5*(X[t2_index][wat+1] + X[t2_index][wat+2]))*0.5;
	    y_avg[t2_index] = (Y[t2_index][wat] + 0.5*(Y[t2_index][wat+1] + Y[t2_index][wat+2]))*0.5;
	    z_avg[t2_index] = (Z[t2_index][wat] + 0.5*(Z[t2_index][wat+1] + Z[t2_index][wat+2]))*0.5;
	    // take only oxygen
	    /* 	      x_avg[t1_index] = X[t1_index][wat]; */
	    /* 	      y_avg[t1_index] = Y[t1_index][wat]; */
	    /* 	      z_avg[t1_index] = Z[t1_index][wat]; */
	    /* 	      x_avg[t2_index] = X[t2_index][wat]; */
	    /* 	      y_avg[t2_index] = Y[t2_index][wat]; */
	    /* 	      z_avg[t2_index] = Z[t2_index][wat]; */

	    displacement.xx = displacement.xy = displacement.yy = 0.0f;
	    displacement.xz = displacement.yz = displacement.zz = 0.0f;

	    if(x_avg[t1_index]>=x_min && x_avg[t1_index]<=x_max)
	      if(y_avg[t1_index]>=y_min && y_avg[t1_index]<=y_max)
		if(z_avg[t1_index]>=z_min && z_avg[t1_index]<=z_max) {
		  displacement.xx = Dab(x_avg[t1_index],x_avg[t2_index],x_avg[t1_index],x_avg[t2_index],timescale);
		  displacement.xy = Dab(x_avg[t1_index],x_avg[t2_index],y_avg[t1_index],y_avg[t2_index],timescale);
		  displacement.yy = Dab(y_avg[t1_index],y_avg[t2_index],y_avg[t1_index],y_avg[t2_index],timescale);
		  displacement.xz = Dab(x_avg[t1_index],x_avg[t2_index],z_avg[t1_index],z_avg[t2_index],timescale);
		  displacement.yz = Dab(y_avg[t1_index],y_avg[t2_index],z_avg[t1_index],z_avg[t2_index],timescale);
		  displacement.zz = Dab(z_avg[t1_index],z_avg[t2_index],z_avg[t1_index],z_avg[t2_index],timescale);

		  x_index = (int) nearbyint((double)((x_avg[t1_index]-x_min)/res)); if (x_index>=nb_elementsX) x_index = nb_elementsX-1;
		  y_index = (int) nearbyint((double)((y_avg[t1_index]-y_min)/res)); if (y_index>=nb_elementsY) y_index = nb_elementsY-1;
		  z_index = (int) nearbyint((double)((z_avg[t1_index]-z_min)/res)); if (z_index>=nb_elementsZ) z_index = nb_elementsZ-1;

		  if (displacement.xx*2*timescale < S->bx*S->bx)
		    if (displacement.yy*2*timescale < S->by*S->by)
		      if (displacement.zz*2*timescale < S->bz*S->bz)
			if (gridDS[x_index][y_index][z_index] >= cut_off) {
			  tensor_count[x_index][y_index][z_index]++;
			  grid[x_index][y_index][z_index] += (displacement.xx*displacement.xx+displacement.yy*displacement.yy+displacement.zz*displacement.zz);
			  gridXX[x_index][y_index][z_index] += displacement.xx;
			  gridXY[x_index][y_index][z_index] += displacement.xy;
			  gridYY[x_index][y_index][z_index] += displacement.yy;
			  gridXZ[x_index][y_index][z_index] += displacement.xz;
			  gridYZ[x_index][y_index][z_index] += displacement.yz;
			  gridZZ[x_index][y_index][z_index] += displacement.zz;
			}
		}
	  }
	  t1_index = t2_index;
	}
      }
    printf("\n");
    close_dcd();
  }


  for(i=0; i<nb_elementsX; i++)
    for(j=0; j<nb_elementsY; j++)
      for(k=0; k<nb_elementsZ; k++)
	{
          if(tensor_count[i][j][k] == 0.0) tensor_count[i][j][k] = 1.0;
	  grid[i][j][k]   /= (6*timescale*(total_frames-5));
	  grid[i][j][k]   *= 1000;
	  gridO[i][j][k]  /= total_frames;
	  gridO[i][j][k]  *= 1000;
	  gridH[i][j][k]  /= total_frames;
	  gridH[i][j][k]  *= 1000;

	  gridXX[i][j][k] /= tensor_count[i][j][k];
	  gridXY[i][j][k] /= tensor_count[i][j][k];
	  gridYY[i][j][k] /= tensor_count[i][j][k];
	  gridXZ[i][j][k] /= tensor_count[i][j][k];
	  gridYZ[i][j][k] /= tensor_count[i][j][k];
	  gridZZ[i][j][k] /= tensor_count[i][j][k];

	}

  EIGEN*** direct = calloc(nb_elementsX,sizeof(EIGEN**));
  for(i=0; i<nb_elementsX; i++) {
    direct[i] = calloc(nb_elementsY,sizeof(EIGEN*));
    for(j=0; j<nb_elementsY; j++) {
      direct[i][j] = calloc(nb_elementsZ,sizeof(EIGEN));
    }
  }
  printf("\nDiagonalizing tensors..\n");
  //computing diffusion direction by finding eigenvalues & eigenvectors
  int best =0;
  float val;
  float th = (float)sqrt(1.5f);
  float ADC,FA;
  double maxADC, maxTCount;
  float term1,term2,term3;

  maxADC=maxTCount=0.0f;
  for(i=0; i<nb_elementsX; i++)
    for(j=0; j<nb_elementsY; j++)
      for(k=0; k<nb_elementsZ; k++)
	{
	  //lower triangle setup of 3x3 matrix
	  AA[0][0] = gridXX[i][j][k];
	  AA[0][1] = gridXY[i][j][k]; AA[1][1] = gridYY[i][j][k];
	  AA[0][2] = gridXZ[i][j][k]; AA[1][2] = gridYZ[i][j][k]; AA[2][2] = gridZZ[i][j][k];

	  dsyev_(&JOBZ, &UPLO, &N, &AA, &LDA, &W, &WORK, &LWORK, &INFO);

	  val = -1.0e-8f; best=0;
	  for(o=0; o<3; o++)
	    if(W[o] > val) {best = o; val = W[o];}
	  direct[i][j][k].b = best;

	  direct[i][j][k].e[0] = W[0]; direct[i][j][k].egv_x[0] = AA[0][0]; direct[i][j][k].egv_y[0] = AA[1][0]; direct[i][j][k].egv_z[0] = AA[2][0];
	  direct[i][j][k].e[1] = W[1]; direct[i][j][k].egv_x[1] = AA[0][1]; direct[i][j][k].egv_y[1] = AA[1][1]; direct[i][j][k].egv_z[1] = AA[2][1];
	  direct[i][j][k].e[2] = W[2]; direct[i][j][k].egv_x[2] = AA[0][2]; direct[i][j][k].egv_y[2] = AA[1][2]; direct[i][j][k].egv_z[2] = AA[2][2];

	  ADC = (W[0] + W[1] + W[2])/3.0f;
          if(ADC>maxADC)maxADC=ADC;
          if(tensor_count[i][j][k]>maxTCount)maxTCount=tensor_count[i][j][k];
	  direct[i][j][k].ADC = ADC;
	  term1 = (W[0]-ADC); term2 = (W[1]-ADC); term3 = (W[2]-ADC);
	  FA = th * sqrt( (term1*term1 + term2*term2 + term3*term3)/ (W[0]*W[0]+W[1]*W[1]+W[2]*W[2]) );
	  direct[i][j][k].FA = FA;
	}

  if (S->tensors=='1') {
    FILE* ten;
    float grid_tmp;

    //    ten = fopen("tensors.sit","wt");
    ten = fopen("tensors.sit","wb");
    //    fprintf(ten,"%i %i %i\n",nb_elementsX,nb_elementsY,nb_elementsZ);
    //    fprintf(ten,"%3.10f %3.10f %3.10f %3.10f\n",x_min+0.5*res,y_min+0.5*res,z_min+0.5*res,res);
    fwrite(&nb_elementsX,sizeof(int),1,ten);
    fwrite(&nb_elementsY,sizeof(int),1,ten);
    fwrite(&nb_elementsZ,sizeof(int),1,ten);
    float min_tmp[4];
    min_tmp[0]=x_min+0.5*res;min_tmp[1]=y_min+0.5*res;min_tmp[2]=z_min+0.5*res;min_tmp[3]=res;
    fwrite(&min_tmp,sizeof(int),4,ten);

    for(k=0; k<nb_elementsZ; k++) {
      for(j=0; j<nb_elementsY; j++) {
	for(i=0; i<nb_elementsX; i++) {
	  grid_tmp=(float)gridXX[i][j][k];
	  if(isnan(grid_tmp))printf("NaN in the grid!\n");
	  //							  fprintf(ten, "%16.8e",grid_tmp);
	  fwrite(&grid_tmp,sizeof(float),1,ten);
	}
      }
    }
    for(k=0; k<nb_elementsZ; k++) {
      for(j=0; j<nb_elementsY; j++) {
      	for(i=0; i<nb_elementsX; i++) {
	  grid_tmp=(float)gridXY[i][j][k];
	  if(isnan(grid_tmp))printf("NaN in the grid!\n");
	  //           		fprintf(ten, "%16.8e  ",grid_tmp);
	  fwrite(&grid_tmp,sizeof(float),1,ten);
	}
      }
    }
    for(k=0; k<nb_elementsZ; k++) {
      for(j=0; j<nb_elementsY; j++) {
      	for(i=0; i<nb_elementsX; i++) {
      	  grid_tmp=(float)gridYY[i][j][k];
	  if(isnan(grid_tmp))printf("NaN in the grid!\n");
	  //       		fprintf(ten, "%16.8e  ",grid_tmp);
	  fwrite(&grid_tmp,sizeof(float),1,ten);
	}
      }
    }
    for(k=0; k<nb_elementsZ; k++) {
      for(j=0; j<nb_elementsY; j++) {
      	for(i=0; i<nb_elementsX; i++) {
      	  grid_tmp=gridXZ[i][j][k];
	  if(isnan(grid_tmp))printf("NaN in the grid!\n");
	  //    		fprintf(ten, "%16.8e  ",grid_tmp);
	  fwrite(&grid_tmp,sizeof(float),1,ten);
	}
      }
    }
    for(k=0; k<nb_elementsZ; k++) {
      for(j=0; j<nb_elementsY; j++) {
      	for(i=0; i<nb_elementsX; i++) {
      	  grid_tmp=(float)gridYZ[i][j][k];
	  if(isnan(grid_tmp))printf("NaN in the grid!\n");
	  //	fprintf(ten, "%16.8e  ",grid_tmp);
	  fwrite(&grid_tmp,sizeof(float),1,ten);
    	}
      }
    }
    for(k=0; k<nb_elementsZ; k++) {
      for(j=0; j<nb_elementsY; j++) {
	for(i=0; i<nb_elementsX; i++) {
	  grid_tmp=gridZZ[i][j][k];
	  if(isnan(grid_tmp))printf("NaN in the grid!\n");
	  //        		fprintf(ten, "%16.8e  ",grid_tmp);
	  fwrite(&grid_tmp,sizeof(float),1,ten);
	}
      }
    }
    fclose(ten);
  }

  printf("\n\n             ** IMPORTANT **\n ADC values in ADC.pdb are scaled (max_ADC = 99.99)  \n Multiply ADC by %e to undo scaling\n", maxADC/99.99);
  if (S->diff=='1') fp1 = fopen("diff.pdb","wt+");
  if (S->oxy=='1') fp2 = fopen("OxygenOUTPDB.pdb","wt+");
  if (S->hydro=='1') fp3 = fopen("HydrogenOUTPDB.pdb","wt+");
  if (S->adc=='1') {
    fp4 = fopen("ADC.pdb","wt+");
    fprintf(fp4,"ADC correction factor = %e\n", maxADC/99.99);
  }
  if (S->fa=='1') fp5 = fopen("FA.pdb","wt+");
  printf("\n");

  int count1=0;
  //SAVING THE CUBE'S CORNER SO WE HAVE SAME OUTPUT - needed for merging of Oxygen.sit and Hydrogen.sit

  //1st corner (bottom)
  x_index = 0; y_index = 0; z_index = 0;
  if (S->diff=='1') {fprintf(fp1,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)1,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->oxy=='1')  {fprintf(fp2,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)1,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->hydro=='1')  {fprintf(fp3,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)1,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->adc=='1') { fprintf(fp4,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)1,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->fa=='1') { fprintf(fp5,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)1,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}

  //2nd corner (bottom)
  x_index = nb_elementsX-1; y_index = 0; z_index = 0;
  if (S->diff=='1') {fprintf(fp1,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)2,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->oxy=='1')  {fprintf(fp2,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)2,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->hydro=='1') {fprintf(fp3,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)2,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->adc=='1') { fprintf(fp4,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)2,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->fa=='1') { fprintf(fp5,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)2,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}

  //3rd corner (bottom)
  x_index = 0; y_index = nb_elementsY-1; z_index = 0;
  if (S->diff=='1') {
    fprintf(fp1,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)3,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->oxy=='1')  {fprintf(fp2,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)3,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->hydro=='1') {fprintf(fp3,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)3,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->adc=='1') { fprintf(fp4,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)3,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->fa=='1') { fprintf(fp5,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)3,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}

  //4th corner (bottom)
  x_index = nb_elementsX-1; y_index = nb_elementsY-1; z_index = 0;
  if (S->diff=='1') {fprintf(fp1,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)4,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->oxy=='1') {fprintf(fp2,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)4,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->hydro=='1')  {fprintf(fp3,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)4,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->adc=='1') { fprintf(fp4,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)4,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->fa=='1') { fprintf(fp5,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)4,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}

  //1st corner (top)
  x_index = 0; y_index = 0; z_index = nb_elementsZ-1;
  if (S->diff=='1') {fprintf(fp1,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)5,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->oxy=='1') {fprintf(fp2,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)5,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->hydro=='1')  {fprintf(fp3,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)5,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->adc=='1') { fprintf(fp4,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)5,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->fa=='1'){ fprintf(fp5,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)5,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}

  //2nd corner (top)
  x_index = nb_elementsX-1; y_index = 0; z_index = nb_elementsZ-1;
  if (S->diff=='1') {fprintf(fp1,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)6,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->oxy=='1') {fprintf(fp2,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)6,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->hydro=='1') {fprintf(fp3,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)6,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->adc=='1') { fprintf(fp4,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)6,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->fa=='1') { fprintf(fp5,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)6,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}

  //3rd corner (top)
  x_index = 0; y_index = nb_elementsY-1; z_index = nb_elementsZ-1;
  if (S->diff=='1') {fprintf(fp1,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)7,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->oxy=='1') {fprintf(fp2,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)7,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->hydro=='1') {fprintf(fp3,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)7,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->adc=='1') { fprintf(fp4,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)7,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->fa=='1') { fprintf(fp5,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)7,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}

  //4th corner (top)
  x_index = nb_elementsX-1; y_index = nb_elementsY-1; z_index = nb_elementsZ-1;
  if (S->diff=='1') {fprintf(fp1,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)8,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->oxy=='1')  {fprintf(fp2,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)8,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->hydro=='1'){fprintf(fp3,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)8,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->adc=='1')  {fprintf(fp4,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)8,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}
  if (S->fa=='1')   {fprintf(fp5,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)8,"DENS","WAT",(long int)1,x_index*res+x_min+res/2,y_index*res+y_min+res/2,z_index*res+z_min+res/2,0.0,0.0f);}

  //  FILE* script;
  //  script = fopen("vmd_script","wt+");
  float d1,d2,d3;
  for(k=0; k<nb_elementsZ; k++) {
    for(j=0; j<nb_elementsY; j++) {
      for(i=0; i<nb_elementsX; i++) {
	if (S->diff=='1') if(grid[i][j][k]>0.0f) {
	    fprintf(fp1,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)count1+1,"DENS","WAT",(long int)1,i*res+x_min+res/2,j*res+y_min+res/2,k*res+z_min+res/2,grid[i][j][k],0.0f);
	  }
 	if (S->oxy=='1') if(gridO[i][j][k]>0.0f) {
	    fprintf(fp2,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)count1+1,"DENS","WAT",(long int)1,i*res+x_min+res/2,j*res+y_min+res/2,k*res+z_min+res/2,gridO[i][j][k],0.0f);
	  }
	if (S->hydro=='1') if(gridH[i][j][k]>0.0f) {
	    fprintf(fp3,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)count1+1,"DENS","WAT",(long int)1,i*res+x_min+res/2,j*res+y_min+res/2,k*res+z_min+res/2,gridH[i][j][k],0.0f);
	  }
	if (direct[i][j][k].e[direct[i][j][k].b]>0.0f) {

	  if (S->adc=='1') fprintf(fp4,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)count1+1,"DENS","WAT",(long int)1,i*res+x_min+res/2,j*res+y_min+res/2,k*res+z_min+res/2,direct[i][j][k].ADC*99.99/maxADC,tensor_count[i][j][k]*99.99/maxTCount);

	  if (S->fa=='1') fprintf(fp5,"ATOM %6li %s %s %5li    %8.3f%8.3f%8.3f%6.2f%6.2f \n",(long int)count1+1,"DENS","WAT",(long int)1,i*res+x_min+res/2,j*res+y_min+res/2,k*res+z_min+res/2,direct[i][j][k].FA,0.0f);
	  d1 = i*res+x_min+res/2; d2 = j*res+y_min+res/2; d3 = k*res+z_min+res/2;
	  //	  fprintf(script,"draw line {%f %f %f} {%f %f %f}\n",d1,d2,d3, d1+(direct[i][j][k].egv_x[direct[i][j][k].b]*res/2.0f),d2+(direct[i][j][k].egv_y[direct[i][j][k].b]*res/2.0f), d3+(direct[i][j][k].egv_z[direct[i][j][k].b]*res/2.0f));
	}
      }
    }
  }
  if (S->diff=='1') fclose(fp1);
  if (S->oxy=='1') fclose(fp2);
  if (S->hydro=='1') fclose(fp3);
  if (S->adc=='1') fclose(fp4);
  if (S->fa=='1') fclose(fp5);
  //  fclose(script);
}

inline static void welcome() {
  printf(">> Water Density & Diffusion v1.1\n");
  printf(">> type 'help' to see a list of commands\n");
  printf(">> type 'print' to see the current values\n");
  printf(">> \n");
}

inline static void print_help() {
  printf(">> Commands:\n");
  printf(">> BOXX, BOXY, BOXZ [float], maximum allowed water displacements (Angstrom)\n");
  printf(">> XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX [float], region of interest (Angstrom)\n");
  printf(">> DENSITY [float], grid density (1/Angstrom)\n");
  printf(">> TSCALE [float], time between trajectory frames (ps), adjusted by FSTEP! \n");
  printf(">> CUTOFF [float], cutoff on occupancy of grid cells with water \n");
  printf(">> FSTEP [integer], trajectory step, 1 - take each frame\n");
  printf(">> ADC, DIFF, FA, HYDRO, OXY [integer 0|1], trigger saving these calculations.\n   the values are saved in the occupancy field of .pdb files\n");
  printf(">> TENSORS [integer 0|1], trigger saving diffusion tensors in .sit file \n");
  printf(">> PRMTOP [string], amber7 parameter file\n");
  printf(">> LOADDCD [sring], .dcd trajectory file  \n");
  printf(">> HELP, PRINT, GO, QUIT\n");
  printf("\n   DIFF, diffusion coefficient D=<ds*ds>/6\n   ADC, apparent diffusion coefficient\n   FA, fractional anisotropy\n   HYDRO, water hydrogen density\n   OXY, water oxygen density\n");
}

inline static void print_setup(Setup* S) {
  int i;
  printf(">> \n");
  printf(">> SETUP PARAMETERS:\n");
  printf(">> \n");
  printf(">> Box X: %3.8f Box Y: %3.8f Box Z: %3.8f\n",S->bx,S->by,S->bz);
  printf(">> Xmin: %3.8f Xmax: %3.8f\n",S->xmin,S->xmax);
  printf(">> Ymin: %3.8f Ymax: %3.8f\n",S->ymin,S->ymax);
  printf(">> Zmin: %3.8f Zmax: %3.8f\n",S->zmin,S->zmax);
  printf(">> Grid Density: %3.8f\n",S->density);
  printf(">> Time between frames: %3.8f\n",S->timescale);
  printf(">> Cut off: %3.8f\n",S->cutoff);
  printf(">> Trajectory step: %3i\n",S->fstep);
  printf(">> Saving ADC calculations to PDB file: "); if (S->adc=='1') printf("YES\n"); else printf("NO\n");
  printf(">> Saving Diffusion calculations to PDB file: "); if (S->diff=='1') printf("YES\n"); else printf("NO\n");
  printf(">> Saving FA calculations to PDB file: "); if (S->fa=='1') printf("YES\n"); else printf("NO\n");
  printf(">> Saving Hydrogen density calculations to PDB file: "); if (S->hydro=='1') printf("YES\n"); else printf("NO\n");
  printf(">> Saving Oxygen density calculations to PDB file: "); if (S->oxy=='1') printf("YES\n"); else printf("NO\n");
  printf(">> Saving Tensors calculations to SIT file: "); if (S->tensors=='1') printf("YES\n"); else printf("NO\n");
  printf(">> Prmtop: %s\n",S->prmtop);
  printf(">> Dcd files[%i]\n",S->files_count); for(i=0; i<S->files_count; i++) {printf(">> [%i]: %s\n",i,S->dcd_files[i]);}
}

inline static void take_command(Setup* S) {
  char buffer[512];
  char command[512];
  char str[512];
  int i=0;

  if (fgets(buffer,500,stdin)==NULL && feof(stdin)) {
    dup2(1,0);//reset to keyboard if input was piped from input file
    fgets(buffer,500,stdin);
  }
  sscanf(buffer,"%s%s",command,str);

  //convert to upper-case
  while(command[i]) {command[i] = (char)toupper(command[i]); i++;}

  if (strncmp(command,"XMIN",4)==0) S->xmin = (float)atof(str);
  else if (!strncmp(command,"BOXX",4)) S->bx = (float)atof(str);
  else if (!strncmp(command,"BOXY",4)) S->by = (float)atof(str);
  else if (!strncmp(command,"BOXZ",4)) S->bz = (float)atof(str);
  else if (!strncmp(command,"XMAX",4)) S->xmax = (float)atof(str);
  else if (!strncmp(command,"YMIN",4)) S->ymin = (float)atof(str);
  else if (!strncmp(command,"YMAX",4)) S->ymax = (float)atof(str);
  else if (!strncmp(command,"ZMIN",4)) S->zmin = (float)atof(str);
  else if (!strncmp(command,"ZMAX",4)) S->zmax = (float)atof(str);
  else if (!strncmp(command,"DENSITY",7)) S->density = (float)atof(str);
  else if (!strncmp(command,"TSCALE",5)) S->timescale = (float)atof(str);
  else if (!strncmp(command,"CUTOFF",6)) S->cutoff = (float)atof(str);
  else if (!strncmp(command,"FSTEP",5)) S->fstep = atoi(str);
  else if (!strncmp(command,"ADC",3))   S->adc = str[0];
  else if (!strncmp(command,"DIFF",4))  S->diff = str[0];
  else if (!strncmp(command,"FA",2))    S->fa = str[0];
  else if (!strncmp(command,"HYDRO",5)) S->hydro = str[0];
  else if (!strncmp(command,"OXY",3))   S->oxy = str[0];
  else if (!strncmp(command,"TENSORS",7)) S->tensors = str[0];
  else if (!strncmp(command,"PRMTOP",6)) {S->prmtop = calloc(strlen(str),sizeof(char)); strcpy(S->prmtop,str);}
  else if (!strncmp(command,"LOADDCD",7)) {
    if(S->files_count==0) S->dcd_files = (char**)calloc(0,sizeof(char));
    S->dcd_files = (char**)realloc(S->dcd_files,sizeof(char*)*(S->files_count+1));
    S->dcd_files[S->files_count] = calloc(strlen(str)+1,sizeof(char));
    strcpy(S->dcd_files[S->files_count],str);
    S->files_count++;
  }
  else if (strncmp(command,"HELP",4)==0) print_help();
  else if (strncmp(command,"PRINT",4)==0) print_setup(S);
  else if (strncmp(command,"GO",2)==0) {
    printf("Calculating water densities...\n");
    WaterDensity(S);
    printf("Done.\n");
    exit(0);
  }
  else if (strncmp(command,"QUIT",4)==0 || strncmp(command,"EXIT",4)==0) {exit(0);}
  memset(buffer,0,256); memset(command,0,256); memset(str,0,256);
}

int main(int argc, char** argv) {
  Setup* s = malloc(sizeof(Setup));
  // Defaults
  s->bx = 20.0;
  s->by = 20.0;
  s->bz = 20.0;
  s->xmin = -20.0; s->xmax = 20.0;
  s->ymin = -20.0; s->ymax = 20.0;
  s->zmin = -20.0; s->zmax = 20.0;
  s->density = 1.0;
  s->timescale = 5.0;
  s->cutoff = 0.001;
  s->fstep = 1.0;
  s->adc =     '1';
  s->diff =    '1';
  s->fa =      '1';
  s->hydro =   '1';
  s->oxy =     '1';
  s->tensors = '1';

  welcome();
  while(1)
    take_command(s);
  return(0);
}
