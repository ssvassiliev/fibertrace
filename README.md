# README #

### Software for analysis of solvent flow in molecular dynamics simulations ###

Streamline analysis is performed in two steps. First 3D diffusion tensor field is calculated using the program "tensors". Next streamline analysis in spirit of MRI fiber tractography is performed using the "streamline" program. 

### Input: ###
1. Molecular parameters in amber7 format and dcd formatted trajectory.
2. User configurable parameters: tfield.conf and streamline.conf

### Output of the "tensors" program: ###
1. 3D map of apparent diffusion coefficient:    "ADC.pdb" 
2. 3D map of fractional anisotropy:             "FA.pdb"  
3. 3D water density map:                        "Oxygen.pdb", "Hydrogen.pdb"
4. 3D diffusion tensor field:                   "tensors.sit"

### Output of the "streamline" program: ###
1. Streamlines color-coded by  anisotropy: streamline_A.mol2  
2. Streamlines color-coded by  diffusion:  streamline_D.mol2
3. Streamlines color-coded by  direction:  streamline_XYZ.mol2

Color code is saved in the charge section of mol2 files.