TENSOR mlsfilter(EIGEN eig, P po, int N,  float sigma, float cutoff, int n) {

  // Half width of this filter at 1/e level is sigma 
  // po     - location in 3D   PEVec.e[0]=PEVec.e[1]=PEVec.e[2];
  // N      - order of polynom
  // sigma  - half width (1/e) of Gau - kernel 
  // cutoff - integration cutoff in units of Gau width (2 = 1.8%)
  // n      - number of integration points (on semi-axis)

  int l,m,zci;
  double x_, y_, z_;

  int i, j, k, ix, iy, iz, zc, Q;
  int C = N+1;
  int N2 = C*C;
  int NC = N2*C; 
  int *nmp, *rst;
  // Integration and weights
  double x, y, z, x_inc, y_inc, z_inc, vol;
  double x_range, y_range, z_range;
  double  G_xyz, wx, MA;
  double *gw, *px, *py, *pz;
  P point; //new point to interpolate tensor at
  TENSOR T; //interpolated tensor
 
   // LAPACK ROUTINE DPOSV 
  char UPLO = 'L';
  int INFO;
  int NRHS = 6;
  int LDA = NC;
  int LDB = NC; 
  double B[6][NC]; 
  double A[NC][LDA];

  zc = 0;
  for(i=0; i<NC; i++) {
    for(j=0; j<NC; j++) A[i][j] = 0.0f;
    for(j=0;j<NRHS;j++) B[j][i] = 0.0f;
  }
  px=calloc(C,sizeof(double)); 
  py=calloc(C,sizeof(double)); 
  pz=calloc(C,sizeof(double));
  nmp=calloc(3,sizeof(int));
  rst=calloc(3,sizeof(int));
  gw=calloc(2*n+1,sizeof(double)); // Gaussian weights 
  // Minimal allowed filter width
  if(eig.e[0] < 1e-4) eig.e[0] = 1e-4;
  if(eig.e[1] < 1e-4) eig.e[1] = 1e-4;
  if(eig.e[2] < 1e-4) eig.e[2] = 1e-4;
 // Integration limits:  half width of the filter at 1/e level
  x_range =  sigma;
  y_range =  sigma*sqrt(eig.e[1]/eig.e[2]);
  z_range =  sigma*sqrt(eig.e[0]/eig.e[2]);
   // Increments: we scale integration volume by cutoff factor 
  x_inc = cutoff*x_range/n;
  y_inc = cutoff*y_range/n;
  z_inc = cutoff*z_range/n;
  vol = x_inc*y_inc*z_inc;
  // Make array of n Gaussian weights
  wx = 1 / (x_range*x_range); 
  int gt = 2*n+1;

  for(ix = 0, x = -n*x_inc; ix < gt; ix++, x += x_inc) gw[ix] = exp(-x*x*wx);

  // INNER LOOP over integration volume
  for(iz = 0, z = -n*z_inc; iz < gt; iz++, z += z_inc) 
    for(iy = 0, y = -n*y_inc; iy < gt; iy++, y += y_inc) 
      for(ix = 0, x = -n*x_inc; ix < gt; ix++, x += x_inc) {
	G_xyz = gw[ix]*gw[iy]*gw[iz]*vol;
	// Current point in world coordinates
	// Translate and rotate (multiply by eigenvec matrix inverse)      
	point.x = po.x + (eig.egv_x[2]*x + eig.egv_x[1]*y + eig.egv_x[0]*z);
	point.y = po.y + (eig.egv_y[2]*x + eig.egv_y[1]*y + eig.egv_y[0]*z);
	point.z = po.z + (eig.egv_z[2]*x + eig.egv_z[1]*y + eig.egv_z[0]*z);
	// Powers of polynom
	px[0]=py[0]=pz[0]=1.0f;
        for(i=1;i<C;i++)
          for(j=1;j<C;j++)
            for(k=1;k<C;k++){
              pz[i]=pz[i-1]*z;
              py[j]=py[j-1]*y;
              px[k]=px[k-1]*x;
            }
	// Accumulate matrices M and B
	T = interpl_tens(point); 
	if (T.flag) {
	  zc++;
	  for(i=0; i<NC; i++) {
	    Q = i; rst[2] = Q/N2; Q -= rst[2]*N2; rst[1] = Q/C; Q -= rst[1]*C;
	    MA = G_xyz*px[Q]*py[rst[1]]*pz[rst[2]];
	    B[0][i] += T.c[0][0]*MA; B[1][i] += T.c[1][0]*MA;
	    B[2][i] += T.c[1][1]*MA; B[3][i] += T.c[2][0]*MA;
	    B[4][i] += T.c[2][1]*MA; B[5][i] += T.c[2][2]*MA;
	    for(j=0; j<=i; j++) {
	      Q = j; nmp[2] = Q/N2; Q -= nmp[2]*N2; nmp[1] = Q/C; Q -= nmp[1]*C;
	      A[j][i] += MA*px[Q]*py[nmp[1]]*pz[nmp[2]];
	    }
	  }
	}
	/* 	if(G_xyz>0.01) */
	/* 	  fprintf(stdout, "ATOM  %5li %s %s %5i    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
		  (long)0, " CA ", "STR",0, x, y, z, 10*G_xyz, 10*G_xyz); */ 
      }
  // END OF INNER LOOP  
  // Solve system of normal equations 
  dposv_(&UPLO, &NC, &NRHS, &A, &LDA, &B, &LDB, &INFO);
  // Get fitted tensor elements at po from zero power coefficients
  T.c[0][0] = B[0][0];
  T.c[1][0] = B[1][0]; T.c[1][1] = B[2][0];
  T.c[2][0] = B[3][0]; T.c[2][1] = B[4][0]; T.c[2][2] = B[5][0];
  if(zc > n*n*n*0.6) T.flag = 1;
  else  T.flag = 0;
  //printf("%f %f %f\n%f %f %f\n%f %f %f\n",T.c[0][0],T.c[0][1],T.c[0][2],T.c[1][0],T.c[1][1],T.c[1][2],T.c[2][0],T.c[2][1],T.c[2][2]);
  free(gw); 
  free(nmp);free(rst);
  free(px); free(py); free(pz);
  return T;
}


inline double xyzmap1000(float *v) {
  double x,y,z,max,sc;
 
  x = fabs(v[0]);  y = fabs(v[1]);  z = fabs(v[2]);
  max = x;  if(y > x) max = y; if(z > max) max = z;
  sc = 9.4/max;
  return((rint(x*sc) + rint(y*sc)*10 + rint(z*sc)*100)*0.001);
}

inline float mountain(EIGEN E) {
	float denom = (E.e[2] + E.e[1] + E.e[0]);
	if (denom<1e-10) return 0.0f;
	return (E.e[2] - E.e[1]) / denom;
}

inline float* direction(EIGEN E) {
        float* dir = calloc(3,sizeof(float));
	// Eigenvalues are sorted in ascending order
        dir[0] = E.egv_x[2];
        dir[1] = E.egv_y[2];
        dir[2] = E.egv_z[2];
        return dir;
}

void inline invert_vector(float *v1){
  v1[0]*=-1; v1[1]*=-1; v1[2]*=-1;
}


inline float dot_product(float *v1, float *v2) {
  return(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);
}



EIGEN eig_lapack(TENSOR T) {
  // Computing eigenvalues & eigenvectors of diffusion tensor
  EIGEN direct;
  char JOBZ = 'V';
  char UPLO = 'L';
  int N = 3;
  double A[3][3];
  int LDA = 3;
  double W[3];
  double WORK[4*3];
  int LWORK = 4*3;
  int INFO;
  int i,j;  
       
  for(i=0; i<3; i++) 
    for(j=0; j<=i; j++) 
      A[j][i] = T.c[i][j];  
  dsyev_(&JOBZ, &UPLO, &N, &A, &LDA, &W, &WORK, &LWORK, &INFO);
  for(i=0; i<3; i++) {
    direct.e[i] = (float)W[i];
    direct.egv_x[i] = (float)A[i][0];
    direct.egv_y[i] = (float)A[i][1];
    direct.egv_z[i] = (float)A[i][2];
  }  
  return direct;
}



TENSOR interpl_tens(P po) {
// Compute tensor at point po  using trilinear interpolation
  TENSOR tmp;
  int i,j,k,l,m,zc;
  double x_, y_, z_;

  // Move point to the origin of the grid
  po.x = po.x - x_min;
  po.y = po.y - y_min;
  po.z = po.z - z_min;
  
  //If  point is outside our tensor field - reset flag and  return
  tmp.flag=0;
  if (( po.x <  0) |  (po.y <  0) | (po.z <  0) |
      (po.x >=  res*(nb_elementsX-1) ) | (po.y >=  res*(nb_elementsY-1) ) | (po.z >=  res*(nb_elementsZ-1) ))
    return tmp;
  
  // Calculate indexes of the voxel  
  i = (int) floor((double)(po.x/res)); 
  j = (int) floor((double)(po.y/res)); 
  k = (int) floor((double)(po.z/res)); 

  // If some of the corners are undefined  - reset flag and  return
  zc=(fields[i][j][k].flag+
       fields[i+1][j][k].flag+
       fields[i][j+1][k].flag+
       fields[i][j][k+1].flag+
       fields[i+1][j+1][k].flag+
       fields[i+1][j][k+1].flag+
       fields[i][j+1][k+1].flag+
       fields[i+1][j+1][k+1].flag);
  if(zc < 4) return tmp;

  // Otherwise proceed:
  // Move point to the origin of the voxel 
  po.x = po.x  - i*res ;
  po.y = po.y  - j*res ;
  po.z = po.z  - k*res;

  // Interpolate tensor 
  x_ = po.x/res; y_ = po.y/res; z_ = po.z/res;
  for(l=0; l<3; l++)
    for(m=0; m<=l; m++)
      tmp.c[l][m]  = fields[i][j][k].c[l][m]*(1-x_)*(1-y_)*(1-z_) +
	fields[i+1][j][k].c[l][m]*po.x*(1-y_)*(1-z_) +
	fields[i][j+1][k].c[l][m]*(1-x_)*po.y*(1-z_) +
	fields[i][j][k+1].c[l][m]*(1-x_)*(1-y_)*po.z +
	fields[i+1][j][k+1].c[l][m]*po.x*(1-y_)*po.z +
	fields[i][j+1][k+1].c[l][m]*(1-x_)*po.y*po.z +
	fields[i+1][j+1][k].c[l][m]*po.x*po.y*(1-z_) +
	fields[i+1][j+1][k+1].c[l][m]*po.x*po.y*po.z;
 tmp.flag = 1;
 return tmp;
}


inline void init_sample(void) {
  // Make model tensor field
  int x,y,z,i,j,k;
  float dx,dy,dz;

  x=y=z=3;
  nb_elementsX = x; nb_elementsY = y; nb_elementsZ = z;

  x_min = y_min = z_min = 0.0f;
  res = 1.0f;
	
  fields = calloc(x,sizeof(TENSOR**));
  gridXX = calloc(x,sizeof(float**));
  gridXY = calloc(x,sizeof(float**)); 
  gridYY = calloc(x,sizeof(float**));
  gridXZ = calloc(x,sizeof(float**)); 
  gridYZ = calloc(x,sizeof(float**)); 
  gridZZ = calloc(x,sizeof(float**));

  float xtmp[3][3][3];
  float ytmp[3][3][3];
  float ztmp[3][3][3];
        
  for(i=0; i<3; i++)
    for(j=0; j<3; j++) {
      xtmp[i][j][0] = 0.0;
      ytmp[i][j][0] = 0.0;
      ztmp[i][j][0] = 1.0;
    }
  for(i=0; i<3; i++) 
    for(j=0; j<3; j++) {
      xtmp[i][j][1] = -0.5;
      ytmp[i][j][1] = 0.9;
      ztmp[i][j][1] = 0.3;
    }
  for(i=0; i<3; i++) 
    for(j=0; j<3; j++) {
      xtmp[i][j][2] = -0.5;
      ytmp[i][j][2] = -0.1;
      ztmp[i][j][2] = 1.0;
    }

  for(i=0; i<x; i++) {
    fields[i] = calloc(y,sizeof(TENSOR*));
    gridXX[i] = calloc(y,sizeof(float*));
    gridXY[i] = calloc(y,sizeof(float*));
    gridYY[i] = calloc(y,sizeof(float*));
    gridXZ[i] = calloc(y,sizeof(float*)); 
    gridYZ[i] = calloc(y,sizeof(float*)); 
    gridZZ[i] = calloc(y,sizeof(float*));
    for(j=0; j<y; j++) {
      fields[i][j] = calloc(z,sizeof(TENSOR));
      gridXX[i][j] = calloc(z,sizeof(float));
      gridXY[i][j] = calloc(z,sizeof(float)); 
      gridYY[i][j] = calloc(z,sizeof(float));
      gridXZ[i][j] = calloc(z,sizeof(float)); 
      gridYZ[i][j] = calloc(z,sizeof(float)); 
      gridZZ[i][j] = calloc(z,sizeof(float));
    }
  }
  int m,n;

  printf("Constructing %ix%ix%i tensor field:\n",x,y,z);
  for(k=0; k<z; k++) 
    for(j=0; j<y; j++) 
      for(i=0; i<x; i++) {

	dx = xtmp[i][j][k];
	dy = ytmp[i][j][k];
	dz = ztmp[i][j][k];

	gridXX[i][j][k] = dx*dx;
	gridYY[i][j][k] = dy*dy;
	gridZZ[i][j][k] = dz*dz;
	gridXY[i][j][k] = dx*dy; 
	gridXZ[i][j][k] = dx*dz; 
	gridYZ[i][j][k] = dy*dz; 

	fields[i][j][k].flag=1;
	fields[i][j][k].c[0][0] = gridXX[i][j][k]; fields[i][j][k].c[0][1] = gridXY[i][j][k]; fields[i][j][k].c[0][2] = gridXZ[i][j][k];
	fields[i][j][k].c[1][0] = gridXY[i][j][k]; fields[i][j][k].c[1][1] = gridYY[i][j][k]; fields[i][j][k].c[1][2] = gridYZ[i][j][k];
	fields[i][j][k].c[2][0] = gridXZ[i][j][k]; fields[i][j][k].c[2][1] = gridYZ[i][j][k]; fields[i][j][k].c[2][2] = gridZZ[i][j][k];
	      
	fprintf(stderr,"Tensor [%i %i %i]:\n",i,j,k); fflush(stderr);
	for(m=0;m<3;m++) {
	  for(n=0;n<3;n++) 
	    fprintf(stderr,"%8.3f ",fields[i][j][k].c[m][n]);
	  fprintf(stderr,"\n");
	}     
      }
  fprintf(stderr,"\n");
}

inline void read_tensor_file(char* file) {
  FILE* ten;
  int x,y,z,i,j,k;

  ten = fopen(file,"rb");
  if(ten==NULL)
    {printf("** E ** tensor file not found\n"); exit(0);}

  fscanf(ten,"%i %i %i\n",&x,&y,&z);
  fscanf(ten,"%f %f %f %f\n",&x_min,&y_min,&z_min,&res);
  nb_elementsX = x; nb_elementsY = y; nb_elementsZ = z;
 	
  fields = calloc(x,sizeof(TENSOR**));
  gridXX = calloc(x,sizeof(float**));
  gridXY = calloc(x,sizeof(float**)); 
  gridYY = calloc(x,sizeof(float**));
  gridXZ = calloc(x,sizeof(float**)); 
  gridYZ = calloc(x,sizeof(float**)); 
  gridZZ = calloc(x,sizeof(float**));
  for(i=0; i<x; i++) {
    fields[i] = calloc(y,sizeof(TENSOR*));
    gridXX[i] = calloc(y,sizeof(float*));
    gridXY[i] = calloc(y,sizeof(float*));
    gridYY[i] = calloc(y,sizeof(float*));
    gridXZ[i] = calloc(y,sizeof(float*)); 
    gridYZ[i] = calloc(y,sizeof(float*)); 
    gridZZ[i] = calloc(y,sizeof(float*));
    for(j=0; j<y; j++) {
      fields[i][j] = calloc(z,sizeof(TENSOR));
      gridXX[i][j] = calloc(z,sizeof(float));
      gridXY[i][j] = calloc(z,sizeof(float));	
      gridYY[i][j] = calloc(z,sizeof(float));
      gridXZ[i][j] = calloc(z,sizeof(float));	
      gridYZ[i][j] = calloc(z,sizeof(float));	
      gridZZ[i][j] = calloc(z,sizeof(float));
    }
  }

  for(k=0; k<z; k++) 
    for(j=0; j<y; j++) 
      for(i=0; i<x; i++) 
	fread(&gridXX[i][j][k],sizeof(float),1,ten);
  for(k=0; k<z; k++) 
    for(j=0; j<y; j++) 
      for(i=0; i<x; i++) 
	fread(&gridXY[i][j][k],sizeof(float),1,ten); 
  for(k=0; k<z; k++) 
    for(j=0; j<y; j++) 
      for(i=0; i<x; i++) 
	fread(&gridYY[i][j][k],sizeof(float),1,ten);
  for(k=0; k<z; k++) 
    for(j=0; j<y; j++) 
      for(i=0; i<x; i++) 
	fread(&gridXZ[i][j][k],sizeof(float),1,ten);
  for(k=0; k<z; k++) 
    for(j=0; j<y; j++) 
      for(i=0; i<x; i++) 
	fread(&gridYZ[i][j][k],sizeof(float),1,ten);
  for(k=0; k<z; k++) 
    for(j=0; j<y; j++) 
      for(i=0; i<x; i++) 
	fread(&gridZZ[i][j][k],sizeof(float),1,ten);

  for(k=0; k<z; k++) {
    for(j=0; j<y; j++) {
      for(i=0; i<x; i++) {
	if (*(int*)&gridXX[i][j][k]==0)fields[i][j][k].flag=0;
	else fields[i][j][k].flag=1;
	fields[i][j][k].c[0][0] = gridXX[i][j][k]; fields[i][j][k].c[0][1] = gridXY[i][j][k]; fields[i][j][k].c[0][2] = gridXZ[i][j][k];
	fields[i][j][k].c[1][0] = gridXY[i][j][k]; fields[i][j][k].c[1][1] = gridYY[i][j][k]; fields[i][j][k].c[1][2] = gridYZ[i][j][k];
	fields[i][j][k].c[2][0] = gridXZ[i][j][k]; fields[i][j][k].c[2][1] = gridYZ[i][j][k]; fields[i][j][k].c[2][2] = gridZZ[i][j][k];
      }
    }
  }
  fclose(ten);
}



void save_interp_grid(int interp) {
  FILE *script, *fp_out; 
  EIGEN g; 
  TENSOR  tt;	
  P po;
  double d1,d2,d3;
  int i,j,k,count;

  count=((nb_elementsX-1)*interp+1)*((nb_elementsY-1)*interp+1)*((nb_elementsZ-1)*interp+1);
  script = fopen("interp_tensor.vmd","wt+");
  fp_out = fopen("interp_grid.gro","wt+");
  fprintf(fp_out,"%s\n%i\n","Interpolated TENSOR grid",count);
  
  for(k=0; k<((nb_elementsZ-1)*interp+1); k++)
    for(j=0; j<((nb_elementsY-1)*interp+1); j++)
      for(i=0; i<((nb_elementsX-1)*interp+1); i++)
	{
	  d1 = (i*res/interp+x_min);
	  d2 = (j*res/interp+y_min);
	  d3 = (k*res/interp+z_min);
	  // The last points
	  if(i==(nb_elementsX-1)*interp) d1-=1e-6;
	  if(j==(nb_elementsY-1)*interp) d2-=1e-6;
	  if(k==(nb_elementsZ-1)*interp) d3-=1e-6;
	  po.x = (float)d1; po.y = (float)d2; po.z = (float)d3;
	  // Write interpolated grid in  GROMACS format
	  fprintf(fp_out,"%5i%5s%5s%5i%10.4f%10.4f%10.4f\n",i*j*k,"ALA","O",i*j*k,d1,d2,d3);
	  tt = interpl_tens(po);
	  if(tt.flag!=0)
	    {
	      g = eig_lapack(tt);
	      // Write principal eigenvector - VMD script// po.x * 10 --- egv_x[2]*3
	      fprintf(script,"draw line {%f %f %f} {%f %f %f}\n",
		      po.x*1 - (g.egv_x[2])*0.4*res/interp,
		      po.y*1 - (g.egv_y[2])*0.4*res/interp,
		      po.z*1 - (g.egv_z[2])*0.4*res/interp,
		      po.x*1 + (g.egv_x[2])*0.4*res/interp,
		      po.y*1 + (g.egv_y[2])*0.4*res/interp,
		      po.z*1 + (g.egv_z[2])*0.4*res/interp);
	    }
	}
  // Finish writing interpolated grid in  GROMACS format
  for(i=0;i<9;i++) fprintf(fp_out,"%10.5f",1.00);
  fprintf(fp_out,"\n");
  fclose(script);
  fclose(fp_out);
}




void save_box(int i_start, int i_end, int j_start, int j_end, int k_start, int k_end )
{
  FILE *fp;

  printf("\nSaving roi.pdb\n");
  fp=fopen("roi.pdb","wt");
  fprintf(fp, "ATOM  %5i %s %s %5i    %8.3f%8.3f%8.3f%6.3f%6.2f \n",
	  1, " CA ", "BOX", 1, i_start*res + x_min, j_start*res + y_min, k_start*res + z_min, 1.0, 0.0);
  fprintf(fp, "ATOM  %5i %s %s %5i    %8.3f%8.3f%8.3f%6.3f%6.2f \n",
	  2, " CA ", "BOX", 1, i_end*res + x_min, j_start*res + y_min, k_start*res + z_min, 1.0, 0.0);
  fprintf(fp, "ATOM  %5i %s %s %5i    %8.3f%8.3f%8.3f%6.3f%6.2f \n",
	  3, " CA ", "BOX", 1, i_start*res + x_min, j_end*res + y_min, k_start*res + z_min, 1.0, 0.0);
  fprintf(fp, "ATOM  %5i %s %s %5i    %8.3f%8.3f%8.3f%6.3f%6.2f \n",
	  4, " CA ", "BOX", 1, i_start*res + x_min, j_start*res + y_min, k_end*res + z_min, 1.0, 0.0);
  fprintf(fp, "ATOM  %5i %s %s %5i    %8.3f%8.3f%8.3f%6.3f%6.2f \n",
	  5, " CA ", "BOX", 1, i_start*res + x_min, j_end*res + y_min, k_end*res + z_min, 1.0, 0.0);
  fprintf(fp, "ATOM  %5i %s %s %5i    %8.3f%8.3f%8.3f%6.3f%6.2f \n",
	  6, " CA ", "BOX", 1, i_end*res + x_min, j_start*res + y_min, k_end*res + z_min, 1.0, 0.0);
  fprintf(fp, "ATOM  %5i %s %s %5i    %8.3f%8.3f%8.3f%6.3f%6.2f \n",
	  7, " CA ", "BOX", 1, i_end*res + x_min, j_end*res + y_min, k_start*res + z_min, 1.0, 0.0);
  fprintf(fp, "ATOM  %5i %s %s %5i    %8.3f%8.3f%8.3f%6.3f%6.2f \n",
	  8, " CA ", "BOX", 1, i_end*res + x_min, j_end*res + y_min, k_end*res + z_min, 1.0, 0.0);
  fprintf(fp, "CONECT%5i%5i\n",1,2);
  fprintf(fp, "CONECT%5i%5i\n",1,3);
  fprintf(fp, "CONECT%5i%5i\n",1,4);
  fprintf(fp, "CONECT%5i%5i\n",2,6);
  fprintf(fp, "CONECT%5i%5i\n",2,7);
  fprintf(fp, "CONECT%5i%5i\n",3,7);
  fprintf(fp, "CONECT%5i%5i\n",3,5);
  fprintf(fp, "CONECT%5i%5i\n",3,1);
  fprintf(fp, "CONECT%5i%5i\n",4,5);
  fprintf(fp, "CONECT%5i%5i\n",4,6);
  fprintf(fp, "CONECT%5i%5i\n",5,8);
  fprintf(fp, "CONECT%5i%5i\n",6,8);
  fprintf(fp, "CONECT%5i%5i\n",7,8);
  fprintf(fp, "END\n");
  fclose(fp);
}


void save_streamline_pdb(P* trace, int *stream_id, long int N) {
  FILE *fp;
  long int i=0;
  long int j,k;
  char label[20];
  int count;
  float maxL, maxA, minL, minA;

  #define maxp 99997

  count=1;

maxL=maxA=0.0;
minL=minA=1000.0;
for(i=0;i<N;i++){
   if(trace[i].a*10>maxL)maxL=trace[i].a*10;
   if(trace[i].c*10>maxA)maxA=trace[i].c*10;
   if(trace[i].a*10<minL)minL=trace[i].a*10;
   if(trace[i].c*10<minA)minA=trace[i].c*10;
}


  sprintf(label,"streamline_%i.pdb",count);
  //printf("Writing file: %s\n",label);
  fp = fopen(label,"wt");  
  fprintf(fp,"%s\n","STREAMLINE");

  for(i=j=0,count=1; i<N; i++,j++)
    {
      if(j>(maxp-1))
	{
	  // Save reference for colormap 
	  fprintf(fp, "ATOM  %5li %s %s %5i    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
		  j, " CA ", "REF", stream_id[i]+1, trace[i].x, trace[i].y, trace[i].z+dt, maxA, maxL);
	  fprintf(fp, "ATOM  %5li %s %s %5i    %8.3f%8.3f%8.3f%6.3f%6.2f \n",
		  j+1, " CA ", "REF", stream_id[i]+2, trace[i].x, trace[i].y, trace[i].z+dt, minA, minL);
	  // Save connectivity
	  for(k=0; k<(maxp-1); k++)
	    if(stream_id[(count-1)*maxp+k]==stream_id[(count-1)*maxp+1+k])
	      fprintf(fp, "CONECT%5li%5li\n",k,k+1);
	  fprintf(fp,"END\n");
	  fclose(fp);
	  count++;
	  sprintf(label,"streamline_%i.pdb",count);
	  //printf("Writing file: %s\n",label);
	  fp = fopen(label,"wt");  
	  fprintf(fp,"%s\n","STREAMLINE");
	  j=0;	  
	}
      fprintf(fp, "ATOM  %5li %s %s %5i    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
	      j, " CA ", "STR",stream_id[i], trace[i].x, trace[i].y, trace[i].z, trace[i].c*10, trace[i].a*10);
    }
  i--;

  // Save reference for colormap 
  fprintf(fp, "ATOM  %5li %s %s %5i    %8.3f%8.3f%8.3f%6.2f%6.2f \n",
	  j, " CA ", "REF", stream_id[i]+1, trace[i].x, trace[i].y, trace[i].z+dt, maxA, maxL);
  fprintf(fp, "ATOM  %5li %s %s %5i    %8.3f%8.3f%8.3f%6.3f%6.2f \n",
	  j+1, " CA ", "REF", stream_id[i]+2, trace[i].x, trace[i].y, trace[i].z+dt, minA, minL);
  // Save connectivity
  for(k=0; k<j; k++)
    if(stream_id[(count-1)*maxp+k]==stream_id[(count-1)*maxp+1+k])
      fprintf(fp, "CONECT%5li%5li\n",k,k+1);
  fprintf(fp,"END\n");
  fclose(fp);
  printf("\n Diffusion color scale: %.2f to %.2f\n", minL, maxL);
  printf("Anisotropy color scale: %.2f to %.2f\n", minA, maxA);
}



void save_streamline_mol2(P* trace, int *stream_id, int stream_count, long int N, char *colorType) {
FILE *fp;
long int i=0;
long int j,k;
char label[20];
int count;
float maxL, maxA, minL, minA;
int *bins;
float binwidth, nbins, hsum, st;

nbins=100; 
bins=calloc(1,nbins*sizeof(int));
count=1;
maxL=maxA=0.0;
minL=minA=1000.0;
for(i=0;i<N;i++){
   if(trace[i].a>maxL)maxL=trace[i].a;
   if(trace[i].a<minL)minL=trace[i].a;
   if(trace[i].b>maxA)maxA=trace[i].b;
   if(trace[i].b<minA)minA=trace[i].b;
}

if(!strcmp(colorType,"diffusion")){
   binwidth=maxL/nbins;
   printf("\n\nDiffusion color scale:\n   min    max    %s", "%\n");
   fp = fopen("streamline_D.mol2","wt");
   for(i=0;i<N;i++){bins[(int)(trace[i].a/binwidth-1)]++;}
   for(i=0, st=0.9, hsum=0.0;i<nbins;i++)
      {
        hsum+=bins[i];
        if(hsum/N>st){printf("%6.2f %6.2f%5.0f%s\n",minL,(i+1)*binwidth,100*hsum/N,"%"); st+=0.02;}
      }
}
if(!strcmp(colorType,"anisotropy")){
   binwidth=maxA/nbins;
   printf("Anisotropy color scale:\n   min    max    %s", "%\n");
   fp = fopen("streamline_A.mol2","wt");
   for(i=0;i<N;i++){bins[(int)(trace[i].b/binwidth-1)]++;}
   for(i=0, st=0.9, hsum=0.0;i<nbins;i++)
      {
        hsum+=bins[i];
        if(hsum/N>st){printf("%6.2f %6.2f%5.0f%s\n",minA,(i+1)*binwidth,100*hsum/N,"%");st+=0.02;}
      }
}

if(!strcmp(colorType,"direction")){
   fp = fopen("streamline_XYZ.mol2","wt");
}




fprintf(fp,"@<TRIPOS>MOLECULE\ngenerated by fibertrace\n");
fprintf(fp,"%li %li %i %i %i\n",N+2 ,N-stream_count, stream_count, 0, 0); 
fprintf(fp,"SMALL\nUSER_CHARGES\n****\nEnergy = 0\n\n");
fprintf(fp,"@<TRIPOS>ATOM\n");

if(!strcmp(colorType,"diffusion")){
  for(i=0,count=1; i<N; i++)
    {
      fprintf(fp, "%8li%5s%10.4f%10.4f%10.4f%5s%8i%5s%12f\n",
	      i+1, "CA", trace[i].x, trace[i].y, trace[i].z, "CA", stream_id[i]+1 ,"STR", trace[i].a);
    }
i--; 
  // Save reference for colormap 
     fprintf(fp, "%8li%5s%10.4f%10.4f%10.4f%5s%8i%5s%12f\n",
     i+2, "CA", trace[i].x, trace[i].y, trace[i].z, "CA", stream_id[i]+2 ,"REF", minL);
      fprintf(fp, "%8li%5s%10.4f%10.4f%10.4f%5s%8i%5s%12f\n",
     i+3, "CA", trace[i].x, trace[i].y, trace[i].z, "CA", stream_id[i]+2 ,"REF", maxL);
}

if(!strcmp(colorType,"anisotropy")){
  for(i=0,count=1; i<N; i++)
    {
      fprintf(fp, "%8li%5s%10.4f%10.4f%10.4f%5s%8i%5s%12f\n",
              i+1, "CA", trace[i].x, trace[i].y, trace[i].z, "CA", stream_id[i]+1 ,"STR", trace[i].b);
    }
i--;
  // Save reference for colormap 
      fprintf(fp, "%8li%5s%10.4f%10.4f%10.4f%5s%8i%5s%12f\n",
     i+2, "CA", trace[i].x, trace[i].y, trace[i].z, "CA", stream_id[i]+2 ,"REF", minA);
      fprintf(fp, "%8li%5s%10.4f%10.4f%10.4f%5s%8i%5s%12f\n",
     i+3, "CA", trace[i].x, trace[i].y, trace[i].z, "CA", stream_id[i]+2 ,"REF", maxA);
}

if(!strcmp(colorType,"direction")){
  for(i=0,count=1; i<N; i++)
    {
      fprintf(fp, "%8li%5s%10.4f%10.4f%10.4f%5s%8i%5s%12f\n",
              i+1, "CA", trace[i].x, trace[i].y, trace[i].z, "CA", stream_id[i]+1 ,"STR", trace[i].c*10);
    }
i--;
  // Save reference for colormap 
      fprintf(fp, "%8li%5s%10.4f%10.4f%10.4f%5s%8i%5s%12f\n",
     i+2, "CA", trace[i].x, trace[i].y, trace[i].z, "CA", stream_id[i]+2 ,"REF", 0.0);
      fprintf(fp, "%8li%5s%10.4f%10.4f%10.4f%5s%8i%5s%12f\n",
     i+3, "CA", trace[i].x, trace[i].y, trace[i].z, "CA", stream_id[i]+2 ,"REF", 10.22);
}

fprintf(fp,"@<TRIPOS>BOND\n");
  // Save connectivity
  for(k=0,i=0; k<N; k++)
    if(stream_id[k]==stream_id[k+1]){
      fprintf(fp, "%10li%10li%10li%3i\n",i+1,k+1,k+2,1);i++;}
 
fprintf(fp,"@TRIPOS>SUBSTRUCTURE\n1 ****        1 TEMP                        0 ****  **** 0 ROOT\n");
fclose(fp);
}

