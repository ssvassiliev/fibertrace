# README #

### Programs for analysis of solvent flow in molecular dynamics simulations ###
![streamlines.png](https://bitbucket.org/repo/qExpaGG/images/3181802118-streamlines.png)

The flow of water in squalene-hopene cyclase computed from 100 ns long molecular dynamics simulation

#### Workflow ####
The analysis is performed in two steps: 
1. 3D diffusion tensor field is calculated using the program **tfield** 
2. Streamline analysis in the spirit of MRI fiber tractography is performed using the **streamline** program. 

#### Dependencies ####
LAPACK, BLAS
 
#### Input files
* Molecular parameters in amber7 format.
* Molecular dynamics trajectory in dcd format.
* TFIELD parameters.
* STREAMLINE parameters.

#### Required TFIELD configuration parameters
The following parameters are required for every tensor field calculation:
* XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX</br>
**Description:** These keywords define ROI. ROI can be either the whole simulation box or only an essential part of it. Cropping ROI can dramatically speed up calculations when the simulation system is very big.</br> 
**Default Values:** MIN = -20.0, MAX = +20.0

* BOXX, BOXY, BOXZ</br>
**Description:** Maximum allowed water displacement. Some water molecules located close to the box boundary may jump from one side of the box to another between consecutive frames. These waters will have huge velocities and they will create strong artifacts seen as straight lines parallel to axes. Setting BOXX, BOXY, BOXZ somewhat smaller than box dimensions will eliminate these artifacts.</br> 
**Default Value:** 20.0

* PRMTOP, LOADDCD</br>
**Description:** Paths to parameter and trajectory files.</br>

* TSCALE</br> 
**Description:** The time interval between frames in the trajectory multiplied by FSTEP.</br> 
**Default value:** 5.0

### Optional TFIELD parameters 
* DENSITY</br>
**Description:** Grid density (1/Angstrom)</br>
**Default value:**  1.0

* CUTOFF</br>
**Description:** Only grid cells with water occupancy higher than CUTOFF will be used for calculation of tensor field.</br> 
**Default value:** 0.001

* FSTEP</br>
**Description:** If the time interval between trajectory frames is small it is possible to increase it for computation of water displacements. For example, if FSTEP is 1 displacement is calculated from frames 1-0, 2-1, 3-2 ... If FSTEP is 2, it is calculated from frames 2-0, 3-1, 4-2 ...</br> 
**Default value:**  1 

Configuration file <b>streamline.conf</b></br>

    
 Run the programs:</br> 
 ./tensors < tfield.conf</br>
 ./streamline < steamline.conf</br>

Or:</br>
 cd test</br>
 ./run_test.sh 
 
### Getting help with the configuration: ###
Run "tensors" program interactively (without input from configuration files) and type "help" at the prompt.

### Output of the "tensors" program: ###
1. 3D map of apparent diffusion coefficient:    "ADC.pdb" 
2. 3D map of fractional anisotropy:             "FA.pdb"  
3. 3D water density map:                        "Oxygen.pdb", "Hydrogen.pdb"
4. 3D diffusion tensor field:                   "tensors.sit"

The values are saved in the "occupancy" field of the ATOM record. In the case of diffusion coefficient output, the weighting factors (grid occupancy) are saved in the "beta" field. 

### Output of the "streamline" program: ###
1. Streamlines color-coded by  anisotropy: streamline_A.mol2  
2. Streamlines color-coded by  diffusion:  streamline_D.mol2
3. Streamlines color-coded by  direction:  streamline_XYZ.mol2

The color code is saved in the charge section of mol2 files. To visualize streamlines colored by direction use shell script "load_streamlines_mol2.sh" (requires VMD).

### References:###
1. C. Gustafsson, S. Vassiliev, C. Kurten, P.-O. Syren, T. Brinck, MD simulations reveal complex water paths in squalene hopene cyclase - tunnel obstructing mutations alter the movements of water in the active site, ACS Omega. submitted (2017).
2. S. Vassiliev, P. Comte, A. Mahboob, D. Bruce, Tracking the Flow of Water through Photosystem II Using Molecular Dynamics and Streamline Tracing, Biochemistry. 49 (2010) 1873–1881. doi:10.1021/bi901900s.
