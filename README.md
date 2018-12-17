
Table of Contents
=================

   * [Fibertrace](#fibertrace)
      * [Programs for analysis of solvent flow in molecular dynamics simulations](#programs-for-analysis-of-solvent-flow-in-molecular-dynamics-simulations)
         * [Workflow](#workflow)
         * [Dependencies](#dependencies)
         * [Input files](#input-files)
         * [Required TFIELD configuration parameters](#required-tfield-configuration-parameters)
         * [Optional TFIELD parameters](#optional-tfield-parameters)
         * [TFIELD output options](#tfield-output-options)
         * [STREAMLINE configuration parameters](#streamline-configuration-parameters)
         * [Getting help with the configuration:](#getting-help-with-the-configuration)
         * [Output of the “streamline” program:](#output-of-the-streamline-program)
         * [References:](#references)

Created by [gh-md-toc](https://github.com/ekalinin/github-markdown-toc)
# Fibertrace
## Programs for analysis of solvent flow in molecular dynamics simulations
![streamlines.png](https://bitbucket.org/repo/qExpaGG/images/3181802118-streamlines.png)

The flow of water in squalene-hopene cyclase computed from 100 ns long molecular dynamics simulation

### Workflow

The analysis is performed in two steps:

1.  3D diffusion tensor field is calculated using the program **tfield**
2.  Streamline analysis in the spirit of MRI fiber tractography is performed using the **streamline** program.

### Dependencies

LAPACK, BLAS

### Input files

-   Molecular parameters in amber7 format.
-   Molecular dynamics trajectory in dcd format.
-   TFIELD parameters.
-   STREAMLINE parameters.

### Required TFIELD configuration parameters

The following parameters are required for every tensor field calculation:

-   XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX  
    **Description:** These keywords define ROI. ROI can be either the whole simulation box or only an essential part of it. Cropping ROI can dramatically speed up calculations when the simulation system is very big.  
    **Default Values:** MIN = -20.0, MAX = +20.0
    
-   BOXX, BOXY, BOXZ  
    **Description:** Maximum allowed water displacement. Some water molecules located close to the box boundary may jump from one side of the box to another between consecutive frames. These waters will have huge velocities and they will create strong artifacts seen as straight lines parallel to axes. Setting BOXX, BOXY, BOXZ somewhat smaller than box dimensions will eliminate these artifacts.  
    **Default Value:** 20.0
    
-   PRMTOP, LOADDCD  
    **Description:** Paths to parameter and trajectory files.  
    **Default value:** none
    
-   TSCALE  
    **Description:** The time interval between frames in the trajectory multiplied by FSTEP in picoseconds.  
    **Default value:** 5.0
    

### Optional TFIELD parameters

-   DENSITY  
    **Description:** Grid density (1/Angstrom)  
    **Default value:** 1.0
    
-   CUTOFF  
    **Description:** Only grid cells with water occupancy higher than CUTOFF will be used for calculation of tensor field.  
    **Default value:** 0.001
    
-   FSTEP  
    **Description:** If the time interval between trajectory frames is small it is possible to increase it for computation of water displacements. For example, if FSTEP is 1 displacement is calculated from frames 1-0, 2-1, 3-2 … If FSTEP is 2, it is calculated from frames 2-0, 3-1, 4-2 …  
    **Default value:** 1
    

### TFIELD output options

-   ADC  
    **Description:** Save apparent diffusion coefficient map in the occupancy field of the file ADC.pdb. The weighting factors (grid occupancy) are saved in the beta field  
    **Default value:** 1 (YES)
    
-   FA  
    **Description:** Save fractional anisotrory map in the occupancy field of the file FA.pdb  
    **Default value:** 1 (YES)
    
-   DIFF  
    **Description:** Save diffusion map in the occupancy field of the file diff.pdb  
    **Default value:** 1 (YES)
    
-   HYDRO  
    **Description:** Save water oxygen density map in the occupancy field of the file Hydrogen.pdb  
    **Default value:** 1 (YES)
    
-   OXY  
    **Description:** Save water hydrogen density map in the occupancy field of the file Oxygen.pdb  
    **Default value:** 1 (YES)
    
-   TENSORS  
    **Description:** Save diffusion tensors in the file tensors.sit  
    **Default value:** 1 (YES)
    

### STREAMLINE configuration parameters

Run the programs:  
./tensors < tfield.conf  
./streamline < steamline.conf

Or:  
cd test  
./run_test.sh

### Getting help with the configuration:

Run “tensors” program interactively (without input from configuration files) and type “help” at the prompt.

### Output of the “streamline” program:

1.  Streamlines color-coded by anisotropy: streamline_A.mol2
2.  Streamlines color-coded by diffusion: streamline_D.mol2
3.  Streamlines color-coded by direction: streamline_XYZ.mol2

The color code is saved in the charge section of mol2 files. To visualize streamlines colored by direction use shell script “load\_streamlines\_mol2.sh” (requires VMD).

### References:

1.  C. Gustafsson, S. Vassiliev, C. Kurten, P.-O. Syren, T. Brinck, MD simulations reveal complex water paths in squalene hopene cyclase - tunnel obstructing mutations alter the movements of water in the active site, ACS Omega. submitted (2017).
2.  S. Vassiliev, P. Comte, A. Mahboob, D. Bruce, Tracking the Flow of Water through Photosystem II Using Molecular Dynamics and Streamline Tracing, Biochemistry. 49 (2010) 1873–1881. doi:10.1021/bi901900s.
