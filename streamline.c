#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "streamline.h"
#include "streamline_lib.c"
// Compilation with Intel MKL:
// loki: icc fibertracing.c -L/opt/intel/mkl/10.0.1.014/lib/em64t -lmkl -lmkl_lapack -openmp  -o fibertrace -xT
// zeus: icc fibertracing.c -L/opt/intel/mkl/10.1.1.019/lib/em64t -lmkl -lmkl_lapack -openmp  -o fibertrace -xT
// icc fibertracing.c -L/opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64  -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -openmp -lpthread -openmp  -o fibertrace -xT -static
static inline P* fibertrace(P po, int* stream_count, float dt, float min_ani_fiber, float min_fiber_length,\
			    float max_fiber_length, float maxturn, float min_lambda, int filter, int npoly, float gcutoff) {
  P pn, ps;
  P* trace, *trace_tmp;
  TENSOR TIntp;
  EIGEN PEVec, PEVec0, EV;
  int i, j, count, count_tmp;
  float *v_new, *v_old, *v_tmp, pmount;
  float cosMAX = cos(M_PI*maxturn/180);

  //filtering variables
  //only need to declare arrays once
  //the mls_filtering function resets all arrays' elements and components to 0.0f everytime it is called
  //filtering variables

  TIntp = interpl_tens(po);
  if (!TIntp.flag) {return NULL;}

  // Calculate velocity and filter shape at start point
  PEVec = eig_lapack(TIntp);
  if(filter)
    { // Do it twice to avoid "jumps" of direction at the first point
      TIntp = mlsfilter(PEVec, po, npoly, sigma, gcutoff, nb_points);
      PEVec = eig_lapack(TIntp);
      TIntp = mlsfilter(PEVec, po, npoly, sigma, gcutoff, nb_points);
    }

  trace_tmp = calloc(1,sizeof(P));
  trace = calloc(1,sizeof(P));

  count_tmp=1;
  trace_tmp[0] = po;
  ps = po;
  min_fiber_length /= dt;
  max_fiber_length /= dt;

  PEVec0 = PEVec = eig_lapack(TIntp);
  v_old = direction(PEVec0);
  trace_tmp[0].c = (float)xyzmap1000(v_old);
  trace_tmp[0].b = mountain(PEVec0);     // color by anisotropy
  trace_tmp[0].a = PEVec0.e[2];          // color by diff coeff
  // Integrate in one direction
  do
    {
      // Euler step
      /*     pn.x = po.x - dt*v_old[0]; */
      /*     pn.y = po.y - dt*v_old[1]; */
      /*     pn.z = po.z - dt*v_old[2]; */

      // RK2  step
      pn.x = po.x - dt*v_old[0]*0.5;
      pn.y = po.y - dt*v_old[1]*0.5;
      pn.z = po.z - dt*v_old[2]*0.5;

      if(filter) {
	if(count_tmp==1) TIntp = mlsfilter(PEVec0, pn, npoly, sigma, gcutoff, nb_points);
	else TIntp = mlsfilter(PEVec, pn, npoly, sigma, gcutoff, nb_points); }
      else  TIntp = interpl_tens(pn);
      if (!TIntp.flag) break;
      EV = eig_lapack(TIntp); v_tmp = direction(EV);
      if(dot_product(v_old,v_tmp) < 0) invert_vector(v_tmp);
      if(dot_product(v_old,v_tmp) < cosMAX)  {free(v_old); free(v_tmp); break;}
      pn.x = po.x - dt*v_tmp[0];
      pn.y = po.y - dt*v_tmp[1];
      pn.z = po.z - dt*v_tmp[2];
      free(v_tmp);

      // Tensor at new position
      if(filter) TIntp = mlsfilter(PEVec, pn, npoly, sigma, gcutoff, nb_points);
      else TIntp = interpl_tens(pn);

      // Check termination conditions
      if (!TIntp.flag) {free(v_old); break;}
      PEVec = eig_lapack(TIntp);
      pmount = mountain(PEVec);
      if((pmount < min_ani_fiber) || (PEVec.e[2] < min_lambda)) {free(v_old); break;}
      v_new = direction(PEVec);
      if(dot_product(v_old,v_new) < 0) invert_vector(v_new);
      if(dot_product(v_old,v_new) < cosMAX)  {free(v_old); free(v_new); break;}
      // Accept new position
      count_tmp++;
      trace_tmp = realloc(trace_tmp,sizeof(P)*(count_tmp));
      trace_tmp[count_tmp-1] = pn;
      trace_tmp[count_tmp-1].c = (float)xyzmap1000(v_old);
      trace_tmp[count_tmp-1].b = pmount; // color by anisotropy
      trace_tmp[count_tmp-1].a = PEVec.e[2];     // color by diff coeff
      po = pn;
      free(v_old);
      v_old = v_new;
    }
  while(count_tmp < max_fiber_length);
  // Now go from start point in opposite direction
  v_old=direction(PEVec0);
  po=ps;
  count=0;
  do
    {
      // Euler step
      /*     pn.x = po.x + dt*v_old[0]; */
      /*     pn.y = po.y + dt*v_old[1]; */
      /*     pn.z = po.z + dt*v_old[2]; */

      // RK2 step
      pn.x = po.x + dt*v_old[0]*0.5;
      pn.y = po.y + dt*v_old[1]*0.5;
      pn.z = po.z + dt*v_old[2]*0.5;
      if(filter){
	if(!count) TIntp = mlsfilter(PEVec0, pn, npoly, sigma, gcutoff, nb_points);
	else TIntp = mlsfilter(PEVec, pn, npoly, sigma, gcutoff, nb_points);}
      else  TIntp = interpl_tens(pn);
      if (!TIntp.flag) break;
      EV = eig_lapack(TIntp);  v_tmp = direction(EV);
      if(dot_product(v_old,v_tmp) < 0) invert_vector(v_tmp);
      if(dot_product(v_old,v_tmp) < cosMAX) {free(v_old); free(v_tmp); break;}
      pn.x = po.x + dt*v_tmp[0];
      pn.y = po.y + dt*v_tmp[1];
      pn.z = po.z + dt*v_tmp[2];
      free(v_tmp);

      // Tensor at new position
      if(filter) {
	if(!count) TIntp = mlsfilter(PEVec0, pn, npoly, sigma, gcutoff, nb_points);
	else TIntp = mlsfilter(PEVec, pn, npoly, sigma, gcutoff, nb_points);}
      else TIntp = interpl_tens(pn);
      // Check termination conditions
      if (!TIntp.flag) {free(v_old); break;}
      PEVec = eig_lapack(TIntp);
      pmount = mountain(PEVec);
      if((pmount < min_ani_fiber) || (PEVec.e[2] < min_lambda)) {free(v_old); break;}
      v_new = direction(PEVec);
      if(dot_product(v_old,v_new) < 0) invert_vector(v_new);
      if(dot_product(v_old,v_new) < cosMAX) {free(v_old); free(v_new); break;}
      // Accept new position
      count++;
      trace = realloc(trace,sizeof(P)*(count+1));
      trace[count-1] = pn;
      trace[count-1].c = (float)xyzmap1000(v_old);
      trace[count-1].b = pmount;
      trace[count-1].a = PEVec.e[2];
      po = pn;
      free(v_old);
      v_old = v_new;
    }
  while(count + count_tmp < max_fiber_length);

  if((count + count_tmp) <  min_fiber_length) {free(trace_tmp); free(trace); return NULL;}
  // Reorder trace
  trace = realloc(trace, sizeof(P)*(count + count_tmp));
  memmove(&trace[count_tmp], trace, sizeof(P)*count);
  for(i = 0,j = 1;j <= count_tmp; i++, j++)
    trace[i] = trace_tmp[count_tmp - j];
  *stream_count = count + count_tmp;
  // Free  memory
  free(trace_tmp);
  return trace;
}


int main(int argc, char** argv) {
  //read_tensor_file("tensors.sit");
  struct  timeb start, end;
  double elapsed_t;

  P *seed_points, po;
  int i,j,k;
  int count, count_total, stream_count;
  int pt_count = 0;
  //int *stream_id;
  int filter, npoly;
  float gcutoff;

  char usage[]="\n"
    " ------------------------------------------------------------------\n"
    " * This is streamline tracing routine with optional mls filtering  *\n"
    " * Normal usage: streamline < values.in                            *\n"
    " * Tensor field read from file: tensors.sit                        *\n"
    " * Streamlines saved in file: streamline.mol2                      *\n"
    " ------------------------------------------------------------------\n"
    "  p     - print current values of parameters\n"
    "  quit  - exit\n"
    "  go    - start tracing\n"
    " ------------------------------------------------------------------\n\n";
  // DEFAULTS
  res = 1.0f; // Grid resolution
  dt = 0.05f; // Stepsize
  filter = 0; // Filter on/off
  sigma = 0.4f; // MLS halfwidth
  nb_points = 10; // MLS datapoints
  npoly = 1; // Order of polynom
  gcutoff = 2.0f;

  printf("%s",usage);

  float min_ani_seed = 0.5; // starting seed points minimal anisotropy
  float min_ani_fiber = 0.1; // minimal anisotropy
  float min_length = 6; // minimal length of a single fiber
  float max_length = 80; // maximal length of a single fiber
  float max_turn = 70; // maximal angular turn from the previous location of a point in a streamline
  float min_lambda = 0.1; // minimal eigenvalue
  float seed_dens = 2; // density of seed points
  int i_start, i_end;
  int j_start, j_end;
  int k_start, k_end;

  read_tensor_file("tensors.sit");
  //init_sample();
  seed_points = calloc(0,sizeof(P));

  i_start = 0; i_end = nb_elementsX - 1;
  j_start = 0; j_end = nb_elementsY - 1;
  k_start = 0; k_end = nb_elementsZ - 1;

  char str[16];
  float val;
  char lbuf[100];

  while(1)
    {
      fprintf(stderr,"\r>");fflush(stderr);
      fgets(lbuf,80,stdin);
      if(feof(stdin)) {dup2(1,0);}
      sscanf(lbuf,"%s %f\n",str,&val);
      if (strcmp(str,"xmin")==0) {i_start = (int)floor((double)(val-x_min)/res);continue;}
      if (strcmp(str,"xmax")==0) {i_end = (int)floor((double)(val-x_min)/res);continue;}
      if (strcmp(str,"ymin")==0) {j_start = (int)floor((double)(val-y_min)/res);continue;}
      if (strcmp(str,"ymax")==0) {j_end = (int)floor((double)(val-y_min)/res);continue;}
      if (strcmp(str,"zmin")==0) {k_start = (int)floor((double)(val-z_min)/res);continue;}
      if (strcmp(str,"zmax")==0) {k_end = (int)floor((double)(val-z_min)/res);continue;}
      if (strcmp(str,"minani")==0) {min_ani_fiber = val;continue;}
      if (strcmp(str,"minlen")==0) {min_length = val;continue;}
      if (strcmp(str,"maxlen")==0) {max_length = val;continue;}
      if (strcmp(str,"maxturn")==0) {max_turn = val;continue;}
      if (strcmp(str,"minlambda")==0) {min_lambda = val;continue;}
      if (strcmp(str,"dt")==0) {dt = val;continue;}
      if (strcmp(str,"minaniseed")==0){min_ani_seed = val;continue;}
      if (strcmp(str,"seeddens")==0){seed_dens = val;continue;}
      if (strcmp(str,"mlsf")==0){filter = val;continue;}
      if (strcmp(str,"fwidth")==0){sigma = val;continue;}
      if (strcmp(str,"fpoly")==0){npoly = val;continue;}
      if (strcmp(str,"fcutoff")==0){gcutoff = val;continue;}
      if (strcmp(str,"fdata")==0){nb_points = val;continue;}
      if (strcmp(str,"savebox")==0)
	{
	  save_box(i_start, i_end, j_start, j_end, k_start, k_end);
	  continue;
	}
      if (strcmp(str,"p")==0){
	printf("               xmin/xmax ymin/ymax zmin/zmax\n");
	printf("ROI start [%10.2f%10.2f%10.2f] (%i,%i,%i)\n",
	       i_start*res+x_min, j_start*res+y_min, k_start*res+z_min,
	       i_start,j_start,k_start);
	printf("ROI end   [%10.2f%10.2f%10.2f] (%i,%i,%i)\n",
	       i_end*res+x_min,j_end*res+y_min,k_end*res+z_min,
	       i_end,j_end,k_end);
	printf("Tensor field density  %.1f\n", res);
	printf("------------------------------------------\n");
	printf("Minimal anisotropy      [minani]      %.2f\n", min_ani_fiber);
	printf("Minimal fiber length    [minlen]      %.1f\n", min_length);
	printf("Maximal fiber length    [maxlen]      %.1f\n", max_length);
	printf("Maximum turn angle      [maxturn]     %.1f\n", max_turn);
	printf("Minimal eigenvalue      [minlambda]   %.2f\n", min_lambda);
	printf("Integration step size   [dt]          %.2f\n", dt);
	printf("Density of seed points  [seeddens]    %.1f\n",seed_dens );
	printf("Minimal seed height of  [minaniseed]  %.2f\n",min_ani_seed );
	if(filter)
	  {
	    printf("MLS Filtering           [mlsf]         YES\n");
	    printf("MLS half width (1/e)    [fwidth]      %.2f\n",sigma);
	    printf("MLS polynom order       [fpoly]       %i\n",npoly);
	    printf("MLS Gau cutoff (HW)     [fcutoff]     %.2f\n",gcutoff);
	    printf("MLS datapoints (1/e)    [fdata]       %i\n",nb_points);
	  }
	else
	  printf("MLS Filtering           [mlsf]         NO\n");
	printf("------------------------------------------\n");
	continue;}
      if (strcmp(str,"go")==0)break;
      if (strcmp(str,"quit")==0)exit(0);
      else {printf("Unknown keyword\n");continue;}
    }


  save_box(i_start, i_end, j_start, j_end, k_start, k_end);

  if (k_start<0 | j_start<0 | i_start<0) {printf("Invalid region of integration. Now exiting...\n"); exit(0);}
  fprintf(stdout,"Finding seed points, may take a while ..\n");fflush(stdout);
  // Find seed points with high anisotropy
  for(k=k_start; k<k_start+(k_end-k_start)*seed_dens; k++)
    for(j=j_start; j<j_start+(j_end-j_start)*seed_dens; j++)
      for(i=i_start; i<i_start+(i_end-i_start)*seed_dens; i++)
	{
	  po.x=i_start*res + (i-i_start)*res/seed_dens + x_min;
	  po.y=j_start*res + (j-j_start)*res/seed_dens + y_min;
	  po.z=k_start*res + (k-k_start)*res/seed_dens + z_min;
	  if (mountain(eig_lapack(interpl_tens(po))) > min_ani_seed) {
	    seed_points = realloc(seed_points,sizeof(P)*(pt_count+1));
	    seed_points[pt_count]=po;
	    pt_count++;
	  }
	}
  printf("Number of seed points: %i\n",pt_count);
  P* tmp;
  P* ftrace = calloc(0,sizeof(P));
  int* stream_id = calloc(0,sizeof(int));
  count_total = stream_count = count = 0;
  ftime(&start);

#pragma omp parallel for private(tmp,count) shared(stream_count, count_total) schedule (dynamic)
  for(i=0; i<pt_count; i++) {
    fprintf(stderr,"Integrating fibertrace[%i] Total length: %i\r", i, count_total); fflush(stderr);
    tmp = fibertrace(seed_points[i], &count, dt, min_ani_fiber, min_length, max_length, \
		     max_turn, min_lambda, filter, npoly, gcutoff);
    if(tmp!=NULL) {
#pragma omp critical
      {
	ftrace = realloc(ftrace,sizeof(P)*(count + count_total));
	stream_id = realloc(stream_id,sizeof(int)*(count + count_total));
	memmove(&ftrace[count_total],tmp,sizeof(P)*count);
	for(j=0;j<count;j++) stream_id[count_total+j]=stream_count;
	count_total += count; stream_count++;
	free(tmp);
      }
    }
  }
#pragma omp barrier
  //  save_streamline_pdb(ftrace, stream_id, count_total);
  save_streamline_mol2(ftrace, stream_id, stream_count, count_total,"diffusion");
  save_streamline_mol2(ftrace, stream_id, stream_count, count_total,"direction");
  save_streamline_mol2(ftrace, stream_id, stream_count, count_total,"anisotropy");

  ftime(&end);
  elapsed_t=(1000.0*(end.time-start.time)+(end.millitm-start.millitm))/1000;

  printf("\nNumber of streamlines: %i\n",stream_count);
  printf("Cumulative number of points for the set of streamlines: %i\n",count_total);

  fprintf(stderr, "\nJob CPU time: %.3f sec\n", elapsed_t);

  // save_interp_grid(4);
  // save_grid();
  // free(stream_id);
  // free(ftrace);
  // free(seed_points);
  return 0;
}
