/* 	C Library created by Pascal Comte, May 2007
	The purposes of this library are:
	1) Read in a binary .dcd coordinate file and write out a text .crd coordinate file
	In this instance, simply call 'write_dcdfile(char* filename)' where filename is the input .dcd file
	The output will have .crd extension
	2) Read xyz coordinates of all the atoms on a frame-by-frame basis
	In this instance, first you must initialize the file with 'initialize_read_xyz(char* filename)'
	Where filename is the input .dcd file
	After the initialization is completed, you must create x, y and z arrays to store the coordinates
	float x[my_H.N]; float y[my_H.N]; float z[my_H.N]; where my_H is a struct created by calling the
	'initialize_read_xyz(char* filename)' function, and where my_H.N is the total number of atoms
	Then simply call 'read_xyz(frame_num, x, y, z)' and it will store the coordinates of all the atoms
	For the frame number 'frame_num'.
	
	An error return an integer value of '0'
	A success return an integer value of '1'
	It is the programmer's choice to do proper error checking in his/her code.
**/

//#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dcd.h"

H my_H;

int read_dcdheader(char* filename) {
	my_H.charmm = my_H.charmm_extrablock = my_H.charmm_4dims = my_H.NSet = my_H.ISTART = 0;
	my_H.NSAVC = my_H.NAMNF = my_H.N = my_H.FREEINDEXES = 0;
	my_H.fp = NULL;
	my_H.DELTA = 0.0;
	my_H.endoffile = -1;
	
	int i;
	char s[4];
	size_t result;

	my_H.fp = fopen(filename,"rb");

	if (my_H.fp==NULL) {
		fprintf(stdout,"Error: DCD File not found.\n");
		fflush(stdout);
		return 0;
	}

	fseek(my_H.fp, 0, SEEK_END);
	my_H.endoffile = ftell(my_H.fp);
	fprintf(stdout,"----------------- DCD header -----------------\n");
	fprintf(stdout,"File Size: %i\n", my_H.endoffile);
	fflush(stdout);
	//rewind(fp);
	fseek(my_H.fp, 4, SEEK_SET);
	
	result = fread(s,4,1,my_H.fp);
	if ( result == (size_t)1 ) {
	  //		fprintf(stdout,"Cord: %s\n",s);
		fflush(stdout);
	}
	else return 0;
	
	result = fread(&my_H.NSet,4,1,my_H.fp);
	if ( result == (size_t)1 ) {
		fprintf(stdout,"NSet: %i ",my_H.NSet);
		fflush(stdout);
	}
	else return 0;
	
	result = fread(&my_H.ISTART,4,1,my_H.fp);
	if ( result == (size_t)1 ) {
		fprintf(stdout,"ISTART: %i ",my_H.ISTART);
		fflush(stdout);
	}
	else return 0;
	
	result = fread(&my_H.NSAVC,4,1,my_H.fp);	
	if ( result == (size_t)1 ) {
		fprintf(stdout,"NSAVC: %i\n",my_H.NSAVC);
		fflush(stdout);
	}
	else return 0;
	
	//Reposition to 40 from beginning; read namnf, number of free atoms
	//rewind(fp);
	fseek(my_H.fp, 40, SEEK_SET);
	result = fread(&my_H.NAMNF,4,1,my_H.fp);
	if ( result == (size_t)1 ) {
		fprintf(stdout,"Number of free atoms: %i\n",my_H.NAMNF);
		fflush(stdout);
	}
	else return 0;

	//Figure out if we are in CHARMm format
	//rewind(fp);
	fseek(my_H.fp,84,SEEK_SET);
	result = fread(&i,4,1,my_H.fp);	
	if ( result == (size_t)1 ) {
		if ( i == 0 ) {
			my_H.charmm=0;
		}
		else {
			my_H.charmm=1;
			my_H.charmm_extrablock=0;
			my_H.charmm_4dims=0;
			//check for extra block
			//rewind(fp);
			fseek(my_H.fp, 48, SEEK_SET);
			result = fread(&i,4,1,my_H.fp);
			if ( result == (size_t)1 ) {
				if ( i == 1 ) {
					my_H.charmm_extrablock=1;
				}
			}
			else return 0;
			//check for 4dims
			//rewind(fp);
			fseek(my_H.fp,52,SEEK_SET);
			result = fread(&i,4,1,my_H.fp);
			if ( result == (size_t)1 ) {
				if ( i == 1 ) {
					my_H.charmm_4dims=1;
				}
			}
			else return 0;
		}
	}
	else return 0;
	fprintf(stdout,"CHARMm: %i CHARMm Extrablock: %i CHARMm 4dims: %i\n", my_H.charmm, my_H.charmm_extrablock, my_H.charmm_4dims);
	fflush(stdout);
	//Read timestep
	fseek(my_H.fp, 44, SEEK_SET);
	if ( my_H.charmm == 0 ) {
		result = fread(&my_H.DELTA,8,1,my_H.fp);
		if ( result != (size_t)1 ) return 0;
	}
	else {
		result = fread(&my_H.DELTA,4,1,my_H.fp);
		if ( result != (size_t)1 ) return 0;
	}
	fprintf(stdout,"Timestep DELTA: %e\n",my_H.DELTA);
	fflush(stdout);
	//Get the size of the next block, and skip it
	//This is the title
	int newsize,numlines;
	fseek(my_H.fp, 92, SEEK_SET);	
	result = fread(&newsize,4,1,my_H.fp);
	if ( result != (size_t)1 ) return 0;
	result = fread(&numlines,4,1,my_H.fp);
	if ( result != (size_t)1 ) return 0;
	
	fseek(my_H.fp, numlines*80, SEEK_CUR);	
	result = fread(&newsize,4,1,my_H.fp);
	if ( result != (size_t)1 ) return 0;
	
	//Read in a 4, then the number of atoms, then another 4
	result = fread(&i,4,1,my_H.fp);
	if ( result != (size_t)1 ) return 0;
	result = fread(&my_H.N,4,1,my_H.fp);
	if ( result != (size_t)1 ) return 0;
	result = fread(&i,4,1,my_H.fp);
        if ( result != (size_t)1 ) return 0;
	fprintf(stdout,"Number of atoms: %i\n", my_H.N);
	fflush(stdout);
	//Stuff with freeindexes.  Just smile and nod.
	int fsize;
	if ( my_H.NAMNF != 0 ) {
		result = fread(&fsize,4,1,my_H.fp);  //should be N-NAMNF*4
		if ( result !=(size_t)1 ) return 0;
		result = fread(&my_H.FREEINDEXES,4,my_H.N - my_H.NAMNF,my_H.fp);
		if ( result !=(size_t)1 ) return 0;
		result = fread(&fsize,4,1,my_H.fp); //should be N-NAMNF*4
		if ( result !=(size_t)1 ) return 0;
	}	
	fprintf(stdout,"----------------------------------------------\n");
	return 1;
};

int write_crdfile(char* filename) {
	int i = read_dcdheader(filename);
	int j;
	if ( i ) {
		int nsets = my_H.NSet;
		int natoms = my_H.N;
		int pos = 0;
		float x[natoms];
		float y[natoms];
		float z[natoms];
		
		for(i=0; i<natoms; i++) {
			x[i] = 0;
			y[i] = 0;
			z[i] = 0;
		}
		FILE* fout;
		char fout_name[256];
		char fout_new[256];
		char* pch;
		pch = strrchr(filename,'/');
		i = (int)(pch-filename+1);
		j = strlen(filename);

		if ( (i < j) && (i>=0) ) {
			strncpy(fout_name, &filename[i],j-i);
		}
		else {
			strcpy(fout_name,filename);
		}

		pch = strrchr(fout_name,'.');
                i = (int)(pch-fout_name+1);
                j = strlen(fout_name);
                if ( (i >0) && (i<j) ) {
                	strncpy(fout_new,fout_name,i-1);
                }
                else {
                        strcpy(fout_new,fout_name);
                }
		strcat(fout_new,".crd");
		fout = fopen(fout_new,"w+");
		fprintf(stdout,"Output filename: %s\n",fout_new);
		fflush(stdout);
		//free(fout_new);
		//free(fout_name);
		fprintf(fout,"rdparm transformed trajectory\n");
		fflush(fout);

		fprintf(stdout,"Now converting...\n");
		fflush(stdout);
		
		for(i=0; i<nsets; i++) {
			pos = ftell(my_H.fp);
			if ( pos == my_H.endoffile ) break;
			if ( read_dcdstep(x,y,z) == 0 ) {
				printf("Error reading atoms\n");
				return 0;
			}
			for(j=0; j<natoms; j++) {
				fprintf(fout,"%4.3f %4.3f %4.3f ",x[j],y[j],z[j]);
				fflush(fout);
			}
			fprintf(stdout,"\b\b\b\b%i\n",(int)((i+1)/(float)nsets*100));
                        fflush(stdout);
		}
		fprintf(stdout,"\n");
                fflush(stdout);
		fclose(fout);
		close_dcd();
	
	}
	else return 0;//if error return error code
	return 1; //return success code
};

int initialize_read_xyz(char* filename) {
	return read_dcdheader(filename);
};

int reset_file() {
	
	int newsize,numlines,i;
	size_t result;

        fseek(my_H.fp, 92, SEEK_SET);
        result = fread(&newsize,4,1,my_H.fp);
        if ( result != (size_t)1 ) return 0;
        result = fread(&numlines,4,1,my_H.fp);
	if ( result != (size_t)1 ) return 0;

        fseek(my_H.fp, numlines*80, SEEK_CUR);
        result = fread(&newsize,4,1,my_H.fp);
        if ( result != (size_t)1 ) return 0;

        //Read in a 4, then the number of atoms, then another 4
        result = fread(&i,4,1,my_H.fp);
        if ( result != (size_t)1 ) return 0;
        result = fread(&my_H.N,4,1,my_H.fp);
        if ( result != (size_t)1) return 0;
        result = fread(&i,4,1,my_H.fp);
        if ( result != (size_t)1 ) return 0;
       
        //Stuff with freeindexes.  Just smile and nod.
        int fsize;
        if ( my_H.NAMNF != 0 ) {
                result = fread(&fsize,4,1,my_H.fp);  //should be N-NAMNF*4
		if ( result !=(size_t)1 ) return 0;
		result = fread(&my_H.FREEINDEXES,4,my_H.N - my_H.NAMNF,my_H.fp);
		if ( result !=(size_t)1 ) return 0;
      		result = fread(&fsize,4,1,my_H.fp); //should be N-NAMNF*4
		if ( result !=(size_t)1 ) return 0;
        }
       	return 1;
};

int read_xyz(int frame_num, float *x, float *y, float*z) {
	
	if ( reset_file() !=1 ) return 0;
        if ( frame_num < 0 || frame_num > my_H.NSet ) return 0;
	int i;
        int nsets = my_H.NSet;
        int natoms = my_H.N;
        int pos = 0;
        //float x[natoms];
        //float y[natoms];
        //float z[natoms];
        for(i=0; i<natoms; i++) {
        	x[i] = 0;
        	y[i] = 0;
        	z[i] = 0;
        }

  	for(i=0; i<nsets; i++) {
        	pos = ftell(my_H.fp);
        	if ( pos == my_H.endoffile ) break;
        	if ( i != frame_num ) {
			if ( forward_dcdstep() == 0 ) {
        			printf("Error reading atoms\n");
        			return 0;
        		}
		}
		if ( i == frame_num ) {
			if ( read_dcdstep(x,y,z) == 0 ) {
                                printf("Error reading atoms\n");
                                return 0;
                        }
			break;
		}
        }
	return 1;
};

int read_dcdstep(float *x, float *y, float *z) {
	int blocksize;
	size_t result;
	if (my_H.charmm & my_H.charmm_extrablock) {
  		result = fread(&blocksize,4,1,my_H.fp);
		if ( result !=(size_t)1 ) return 0;
		fseek(my_H.fp,blocksize,SEEK_CUR);
  		result = fread(&blocksize,4,1,my_H.fp);
		if ( result !=(size_t)1 ) return 0;
	}	

	if (my_H.NAMNF == 0) {
		//getting x-coordinates
		result = fread(&blocksize,4,1,my_H.fp);
		if ( result !=(size_t)1 ) return 0;
		result = fread(x,4,blocksize/4,my_H.fp);
		if ( result != (size_t)blocksize/4 ) return 0;
       		result = fread(&blocksize,4,1,my_H.fp);
		if ( result !=(size_t)1 ) return 0;
		
		//getting y-coordinates
		result = fread(&blocksize,4,1,my_H.fp);
		if ( result != (size_t)1 ) return 0;
		result = fread(y,4,blocksize/4,my_H.fp);
		if ( result != (size_t)blocksize/4 ) return 0;
       		result = fread(&blocksize,4,1,my_H.fp);
		if ( result !=(size_t)1 ) return 0;
		
		//getting z-coordinates
		result = fread(&blocksize,4,1,my_H.fp);
		if ( result !=(size_t)1 ) return 0;
		result = fread(z,4,blocksize/4,my_H.fp);
		if ( result != (size_t)blocksize/4 ) return 0;
       		result = fread(&blocksize,4,1,my_H.fp);
		if ( result !=(size_t)1 ) return 0;
	}
	else {
		 //this is not implemented in the VMD code I copied from
	}
 
	//Skip the 4th dimension, if present
	if (my_H.charmm & my_H.charmm_4dims) {
 		result = fread(&blocksize,4,1,my_H.fp);
		if ( result !=(size_t)1 ) return 0;
		fseek(my_H.fp, blocksize, SEEK_CUR);
		result = fread(&blocksize,4,1,my_H.fp);
		if ( result !=(size_t)1 ) return 0;
	}			

	return 1;
};

int forward_dcdstep() {
	int blocksize;
	size_t result;
	if (my_H.charmm & my_H.charmm_extrablock) {
  		result = fread(&blocksize,4,1,my_H.fp);
		if ( result !=(size_t)1 ) return 0;
		fseek(my_H.fp,blocksize,SEEK_CUR);
  		result = fread(&blocksize,4,1,my_H.fp);
		if ( result !=(size_t)1 ) return 0;
	}

	if (my_H.NAMNF == 0) {
		//getting x-coordinates
		result = fread(&blocksize,4,1,my_H.fp);
		if ( result !=(size_t)1 ) return 0;
		fseek(my_H.fp,blocksize,SEEK_CUR);
		//result = fread(x,4,blocksize/4,my_H.fp);
		result = fread(&blocksize,4,1,my_H.fp);
		if ( result !=(size_t)1 ) return 0;
		
		//getting y-coordinates
		result = fread(&blocksize,4,1,my_H.fp);
		if ( result !=(size_t)1 ) return 0;
		fseek(my_H.fp,blocksize,SEEK_CUR);
		//result = fread(y,4,blocksize/4,my_H.fp);
		result = fread(&blocksize,4,1,my_H.fp);
		if ( result !=(size_t)1 ) return 0;
			
		//getting z-coordinates
		result = fread(&blocksize,4,1,my_H.fp);
		if ( result !=(size_t)1 ) return 0;
		fseek(my_H.fp,blocksize,SEEK_CUR);
		//result = fread(z,4,blocksize/4,my_H.fp);
		result = fread(&blocksize,4,1,my_H.fp);
		if ( result !=(size_t)1 ) return 0;
	}
	else {
		 //this is not implemented in the VMD code I copied from
	}
 
	//Skip the 4th dimension, if present
	if (my_H.charmm & my_H.charmm_4dims) {
 		result = fread(&blocksize,4,1,my_H.fp);
		if ( result !=(size_t)1 ) return 0;
		fseek(my_H.fp, blocksize, SEEK_CUR);
		result = fread(&blocksize,4,1,my_H.fp);
		if ( result !=(size_t)1 ) return 0;
	}			
	return 1;
};

int close_dcd() {
	if (my_H.fp != NULL) {
		fclose(my_H.fp);
		return 1;
	}
	else {
		printf("Error closing DCD file.\n");
		return 0;
	}
}


