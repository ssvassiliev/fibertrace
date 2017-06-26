# README #

### Programs for analysis of solvent flow in molecular dynamics simulations ###
![streamlines.png](https://bitbucket.org/repo/qExpaGG/images/3181802118-streamlines.png)

Flow of water in squalene hopene cyclase computed from 100 ns long molecular dynamics simulation

### Workflow: ###
Analysis is performed in two steps: first 3D diffusion tensor field is calculated using the program "tensors", next streamline analysis in spirit of MRI fiber tractography is performed using the "streamline" program. 

### Dependencies: ###
LAPACK, BLAS

### Input: ###
1. Molecular parameters in amber7 format and trajectory in dcd format.
2. User configurable parameters: "tfield.conf" and "streamline.conf".

### Running the programs: ###
Edit "tfield.conf" to define your ROI and point to MD simulation files, then run the programs: 

 ./tensors < tfield.conf

 ./streamline < steamline.conf

Or:
 
 cd test

 ./run_test.sh 
 

### Getting help with the configuration: ###
Run "tensors" program interactively (without input from configuration files) and type "help" at the prompt.

### Output of the "tensors" program: ###
1. 3D map of apparent diffusion coefficient:    "ADC.pdb" 
2. 3D map of fractional anisotropy:             "FA.pdb"  
3. 3D water density map:                        "Oxygen.pdb", "Hydrogen.pdb"
4. 3D diffusion tensor field:                   "tensors.sit"

The values are saved in the "occupancy" field of the ATOM record. In case of diffusion coefficient output the weighting factors (grid occupancy) are saved in the "beta" field. 


### Output of the "streamline" program: ###
1. Streamlines color-coded by  anisotropy: streamline_A.mol2  
2. Streamlines color-coded by  diffusion:  streamline_D.mol2
3. Streamlines color-coded by  direction:  streamline_XYZ.mol2

The color code is saved in the charge section of mol2 files. To visualize streamlines colored by direction use shell script "load_streamlines_mol2.sh" (requires VMD).

### Reference:###
T. Zaraiskaya, S. Vassiliev, D. Bruce, Discovering oxygen channel topology in photosystem II using implicit ligand sampling and wavefront propagation. Journal of Computational Science 5 (3), 549-555, 2014, DOI: 10.1016/j.jocs.2013.10.001