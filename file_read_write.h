#define MAX_NATOMS 900000
#define MAX_NRES   400000
#define MAX_FACETS 10000

#include <math.h>
float roundf(float);
double nearbyint(double);

/*---- Include this in compilation: -----*/
//  "/thunderhome/svassili/src/Lib/file_read_write.h"
//  "/thunderhome/svassili/src/Lib/file_read_write.c"

/*----------------- Common -------------------*/
long good_list[MAX_NRES];
int connect_list[MAX_NATOMS][6];
int connect_num[MAX_NATOMS];
/*---------------- PDB  RECORD(S) ---------------*/
long  N;
char   atom[MAX_NATOMS][7];
long   serial[MAX_NATOMS];
char   atom_name[MAX_NATOMS][6];
char   altLoc[MAX_NATOMS][2];
char   resName[MAX_NATOMS][6];
char   chainID[MAX_NATOMS][3];
long   resSeq[MAX_NATOMS],resSeqN[MAX_NATOMS];
char   iCode[MAX_NATOMS][4];
double X[MAX_NATOMS];
double Y[MAX_NATOMS];
double Z[MAX_NATOMS];
float  occupancy[MAX_NATOMS],tempFactor[MAX_NATOMS];
/* PDB RECORD 2 */
long  N2;
char   atom2[MAX_NATOMS][7];
long   serial2[MAX_NATOMS];
char   atom_name2[MAX_NATOMS][6];
char   altLoc2[MAX_NATOMS][2];
char   resName2[MAX_NATOMS][6];
char   chainID2[MAX_NATOMS][3];
long   resSeq2[MAX_NATOMS],resSeqN2[MAX_NATOMS];
char   iCode2[MAX_NATOMS][4];
double X2[MAX_NATOMS];
double Y2[MAX_NATOMS];
double Z2[MAX_NATOMS];
float  occupancy2[MAX_NATOMS],tempFactor2[MAX_NATOMS];
/*-----------------AMBER PRMTOP ---------------*/
long NRES, NRES2, RESIDUE_POINTER[MAX_NRES];
char RESIDUE_LABEL[MAX_NRES][5];
float charge[MAX_NATOMS], mass[MAX_NATOMS];
char Type[MAX_NATOMS][5], Element[MAX_NATOMS][3];
float radius[MAX_NATOMS];

float vdw_rad[MAX_NATOMS], vdw_eps[MAX_NATOMS];
long NBONH, IBH[MAX_NATOMS], JBH[MAX_NATOMS];
long NBONA, IB[MAX_NATOMS], JB[MAX_NATOMS];

long f1[MAX_FACETS],f2[MAX_FACETS],f3[MAX_FACETS];

               
/*------------------ READING FUNCTIONS -------------------*/
void read_namd_binary(char*, unsigned long, double* ,double* ,double*);
void read_amber_coor(char*, long, double*, double*, double*);
void read_pdb(char *);
void read_pdb2(char *);
int read_parm7(char*);
void read_msms_vert(char*, double*, double*, double*);
// Read QM route from file, store it and return pointer  
char *read_route(char *);
/*---------------------------------------------------------*/


/*------------------ WRITING FUNCTIONS --------------------*/
// WRITE NAMD BINARY RESTART FILE
void write_namd_binary(char*, long, double*, double*, double*);
// WRITE INPUT FOR MSMS
void write_xyzr(char *);
// WRITE MCCE-FORMATTED PARAMETER FILES 
void write_mcce(char *);
void write_mcce_ff(char *);
void write_xyzrc(char *);
//  WRITE INPUT FOR JAGUAR
void write_jaguar(int, char *, char *);
//  WRITE INPUT FOR GAUSSIAN
void write_gau(int, char *);

// WRITE PDB
/*-------------- WRITE SELECTED RESIDUES  ---------------*/
// Selected residue numbers are in array good_list[MAX_NRES]
// At_num > 99999 and resnum > 9999 are written in HEX
// Works with the size of up to 65535 residies and/or up to 1048575 atoms
void write_one_pdb_sel(char *); 
// Split input in chinks of 9999 residues,  
// Anum > 99999 are written in HEX
void write_split_pdb_sel(char *);
/*--------------------------------------------------------*/
/*-------------- WRITE ALL RESIDUES  ---------------*/
void write_one_pdb(char *); 
/*------------- NAMD RESTART  MANIPULATION --------------*/  
void insert_wat_namdbin(char *, double*, double*, double*);
void add_wat(double,double,double);

/*------------------- CONVERSION ------------------------*/
void convert_amb_types_to_radii(void);
/* Uses a Table of msms radii */
void wrap_names(void); // Wrap 4-character at_names
void mass_to_element(void); // Convert mass to atomic symbol

/*------------------ CHECK BONDS --------------------*/
/* Finds and prints out long bonds */
void check_bonds(void);

/*------------------ CENTER THE SYSTEM -----------------*/
void center(void);
/* Translates system to bring its center of gravity to the origin */

/*-- REMOVE WATER FROM SIDES AND MIDDLE OF THE MEMBRANE --*/ 
/* start residue number, membrane zcen, membrane z_width, margin */
void remove_bad_wat(long, float, float, float);
