
Table of Contents

   * [Fibertrace: programs for analysis of solvent flow in molecular dynamics simulations](#fibertrace-programs-for-analysis-of-solvent-flow-in-molecular-dynamics-simulations)
      * [Introduction](#introduction)
      * [Workflow](#workflow)
      * [Dependencies](#dependencies)
      * [Input files](#input-files)
      * [Program parameters and output options](#program-parameters-and-output-options)
         * [Required TFIELD configuration parameters](#required-tfield-configuration-parameters)
         * [Optional TFIELD parameters](#optional-tfield-parameters)
         * [TFIELD output options](#tfield-output-options)
         * [STREAMLINE configuration parameters](#streamline-configuration-parameters)
         * [Running the programs:](#running-the-programs)
      * [Getting help with the configuration:](#getting-help-with-the-configuration)
      * [Output of the “streamline” program:](#output-of-the-streamline-program)
      * [References:](#references)

Created by [gh-md-toc](https://github.com/ekalinin/github-markdown-toc)
# Fibertrace: programs for analysis of solvent flow in molecular dynamics simulations

<p align="center">
  <img src="https://bitbucket.org/repo/qExpaGG/images/3181802118-streamlines.png"><br>
  The flow of water in squalene-hopene cyclase computed from 100 ns long molecular dynamics simulation
</p>

## Introduction

Streamline analysis is based on the determination of the 3D velocity field and the subsequent calculation of the motion of the massless fluid element in it. To determine velocity field the simulation system volume is divided into *n* small cubic elements (voxels) with a volume of ~1&#8491;<sup>3</sup> each. For all voxels containing a water molecule, the diffusion tensor elements are calculated according to the Einstein relation:<br><br>

```math
T^{\alpha\beta}=\frac{\langle \alpha (t + \Delta t)- \alpha(t) \rangle \langle \beta(t + \Delta t) - \beta(t) \rangle}{2\Delta t}\qquad \alpha ,\beta = \{x,y,z\}
```

To obtain all nine tensor elements, &#945; and &#946; are sequentially substituted with *x*, *y*, and *z* coordinates of a water molecule. Tensor elements are averaged over the time window of the molecular dynamics run, then eigenvalues  &#955; and eigenvectors &#957; are obtained by diagonalization of the diffusion tensors.
From the eigenvalues apparent diffusion coefficient (*ADC*) and fractional anisotropy (*FA*) are calculated:<br>

```math
ADC=\frac{\lambda_1 +\lambda_2 + \lambda_3}{3}
```
```math
FA=\sqrt{\frac{3}{2}}\sqrt{\frac{\left(\lambda_1-\frac{\lambda_1 + \lambda_2 + \lambda_3}{3}\right)^2 + \left(\lambda_2-\frac{\lambda_1 + \lambda_2 + \lambda_3}{3}\right)^2 + \left (\lambda_3-\frac{\lambda_1 + \lambda_2 + \lambda_3}{3}\right)^2}{\lambda_1^2 + \lambda_2^2 + \lambda_3^2}}
```

*FA* has values between 0 and 1. A value of 0 corresponds to completely isotropic diffusion in all directions, a value of 1 corresponds to diffusion only along one axis.

## Workflow
The analysis is performed in two steps:
1.  3D diffusion tensor field is calculated and diagonalized using the program **tfield**
2.  Streamline analysis in the spirit of MRI fiber tractography is performed using the **streamline** program.

## Dependencies
LAPACK, BLAS

## Input files
-   Molecular parameters in amber7 format.
-   Molecular dynamics trajectory in dcd format.
-   TFIELD parameters.
-   STREAMLINE parameters.

## Program parameters and output options
### Required TFIELD configuration parameters
The following parameters are required for every tensor field calculation:
-   XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX\
    **Description:** These keywords define ROI. ROI can be either the whole simulation box or only an essential part of it. Cropping ROI can dramatically speed up calculations when the simulation system is very big.\
    **Default Values:** MIN = -20.0, MAX = +20.0

-   BOXX, BOXY, BOXZ\
    **Description:** Maximum allowed water displacement. Some water molecules located close to the box boundary may jump from one side of the box to another between consecutive frames. These waters will have huge velocities and they will create strong artifacts seen as straight lines parallel to axes. Setting BOXX, BOXY, BOXZ somewhat smaller than box dimensions will eliminate these artifacts.\
    **Default Value:** 20.0

-   PRMTOP, LOADDCD\
    **Description:** Paths to parameter and trajectory files.\
    **Default value:** none\
-   TSCALE\
    **Description:** The time interval between frames in the trajectory multiplied by FSTEP, picoseconds.\
    **Default value:** 5.0

### Optional TFIELD parameters
-   DENSITY\
    **Description:** Grid density.\
    **Default value:** 1.0 &#8491;<sup>-1</sup>
-   CUTOFF\
    **Description:** Only grid cells with water occupancy higher than CUTOFF will be used for calculation of tensor field.\
    **Default value:** 0.001

-   FSTEP\
    **Description:** If the time interval between trajectory frames is small it is possible to increase it for computation of water displacements. For example, if FSTEP is 1 displacement is calculated from frames 1-0, 2-1, 3-2 … If FSTEP is 2, it is calculated from frames 2-0, 3-1, 4-2 …\
    **Default value:** 1

### TFIELD output options
-   ADC\
    **Description:** Save apparent diffusion coefficient map. Values are saved in the occupancy field of the file ADC.pdb. The weighting factors (grid occupancy) are saved in the beta field\
    **Default value:** 1 (YES)

-   FA\
    **Description:** Save fractional anisotropy map. Values are saved in the occupancy field of the file FA.pdb\
    **Default value:** 1 (YES)

-   DIFF\
    **Description:** Save diffusion coefficient map. Values are saved in the occupancy field of the file diff.pdb. This is the same coefficient as saved in ADC.pdb, but obtained directly from the mean square displacement:

    ```math
    D=\left\langle\frac{\Delta x^{2}+\Delta y^{2}+\Delta z^{2}}{6\Delta t}\right\rangle
    ```
    **Default value:** 1 (YES)


-   HYDRO\
    **Description:** Save water oxygen density map. Values are saved in the occupancy field of the file Hydrogen.pdb\
    **Default value:** 1 (YES)

-   OXY\
    **Description:** Save water hydrogen density map. Values are saved in the occupancy field of the file Oxygen.pdb\
    **Default value:** 1 (YES)\
-   TENSORS\
    **Description:** Save diffusion tensors. Tensors are saved in the file tensors.sit\
    **Default value:** 1 (YES)

### STREAMLINE configuration parameters
-   xmin, xmax, ymin, ymax, zmin, zmax\
    **Description:** These keywords define ROI.\
    **Default Values:** None

-   res\
    **Description:** Tensor field density.\
    **Default value:** 1.0 &#8491;<sup>-1</sup>

-   dt\
    **Description:** Integration step size\
    **Default value:** 0.05

-   minlen, maxlen\
        **Description:** Minimal and maximal allowed length of a streamline. Streamlines shorter than minlen will be discarded. A streamline is terminated when its length reaches maxlen.\
        **Default values:** minlen = 6.0 &#8491;, maxlen = 80.0 &#8491;

-   maxturn\
        **Description:** Maximal allowed angular turn from the previous location of a point in a streamline.\
        **Default value:** 70.0 degrees

-   minaniseed\
        **Description:** Minimal linear anisotropy *c* of seed points.
        ```math
        c=\frac{\lambda_1-\lambda_2}{\lambda_1+\lambda_2+\lambda_3}
        ```

         Seed points are points with *c* > minaniseed. These points are used as a starting points for propagation of streamlines. Lower value of minaniseed will result in more streamlines.\
        **Default value:** 0.5

-   seeddens\
        **Description:** Density of seed points. Higher value of seed_dens will result in more streamlines.\
        **Default value:** 2.0 &#8491;<sup>-1</sup>

-   mlsf\
        **Description:** Turn moving least-squares filtering on/off. Moving least-squares filter approximates tensor locally with a low-degree polynomial to the location, orientation, and history of a streamline. Set mlsf to 1 to enable filtering.\
        **Default value:** 0 (no filter)</sup>

-   fpoly\
        **Description:** Order of polynomial.\
        **Default value:** 1

-   fwidth\
        **Description:** Half width of the Gaussian filter kernel at e<sup>-1</sup> level.\
        **Default value:** 0.4

-   fdata\
        **Description:** The number of integration points on a semi axis.\
        **Default value:** 10

-   fcutoff\
        **Description:** Integration cutoff in units of Gaussian filter kernel width.\
        **Default value:** 2.0

-   savebox\
        **Description:**  Save ROI in the roi.pdb file. This option has no values. \
        **Default:** do not save ROI


### Running the programs:
`./tensors < tfield.conf`\
`./streamline < streamline.conf`

Or:\
`cd test`\
`./run_test.sh`a

## Getting help with the configuration:

Run “tensors” program interactively (without input from configuration files) and type “help” at the prompt.

## Output of the “streamline” program:
Streamlines are saved in mol2 file format as atoms with atom name CA and residue name STR. Each streamline is represented as one residue i.e. same residue ID is assigned to all atoms of a streamline. Additional information is saved in the charge column of mol2 files (the last column). The streamline program saves calculation results in 3 mol2 files:

1.  Fractional anisotropy is saved in the streamline_A.mol2
2.  The largest eigenvalue of the diffusion tensor ( &#955;<sub>1</sub> ) is saved in the streamline_D.mol2.  It shows diffusivity (the rate of diffusion) in the direction tangent to streamlines.
3.  The direction of streamlines is saved in streamline_XYZ.mol2. Unit vector tangent to a streamline at each point is encoded into a floating point number (colour code). Streamlines coloured by direction can be visualized in VMD using shell script “load\_streamlines\_mol2.sh”.

## References:

1.  C. Gustafsson, S. Vassiliev, C. Kurten, P.-O. Syren, T. Brinck, MD simulations reveal complex water paths in squalene-hopene cyclase - tunnel obstructing mutations alter the movements of water in the active site, ACS Omega. submitted (2017).
2.  S. Vassiliev, P. Comte, A. Mahboob, D. Bruce, Tracking the Flow of Water through Photosystem II Using Molecular Dynamics and Streamline Tracing, Biochemistry. 49 (2010) 1873–1881. doi:10.1021/bi901900s.
