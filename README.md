---


---

<h1 id="readme">README</h1>
<h3 id="programs-for-analysis-of-solvent-flow-in-molecular-dynamics-simulations">Programs for analysis of solvent flow in molecular dynamics simulations</h3>
<p><img src="https://bitbucket.org/repo/qExpaGG/images/3181802118-streamlines.png" alt="streamlines.png"></p>
<p>The flow of water in squalene-hopene cyclase computed from 100 ns long molecular dynamics simulation</p>
<h4 id="workflow">Workflow</h4>
<p>The analysis is performed in two steps:</p>
<ol>
<li>3D diffusion tensor field is calculated using the program <strong>tfield</strong></li>
<li>Streamline analysis in the spirit of MRI fiber tractography is performed using the <strong>streamline</strong> program.</li>
</ol>
<h4 id="dependencies">Dependencies</h4>
<p>LAPACK, BLAS</p>
<h4 id="input-files">Input files</h4>
<ul>
<li>Molecular parameters in amber7 format.</li>
<li>Molecular dynamics trajectory in dcd format.</li>
<li>TFIELD parameters.</li>
<li>STREAMLINE parameters.</li>
</ul>
<h3 id="required-tfield-configuration-parameters">Required TFIELD configuration parameters</h3>
<p>The following parameters are required for every tensor field calculation:</p>
<ul>
<li>
<p>XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX<br>
<strong>Description:</strong> These keywords define ROI. ROI can be either the whole simulation box or only an essential part of it. Cropping ROI can dramatically speed up calculations when the simulation system is very big.<br>
<strong>Default Values:</strong> MIN = -20.0, MAX = +20.0</p>
</li>
<li>
<p>BOXX, BOXY, BOXZ<br>
<strong>Description:</strong> Maximum allowed water displacement. Some water molecules located close to the box boundary may jump from one side of the box to another between consecutive frames. These waters will have huge velocities and they will create strong artifacts seen as straight lines parallel to axes. Setting BOXX, BOXY, BOXZ somewhat smaller than box dimensions will eliminate these artifacts.<br>
<strong>Default Value:</strong> 20.0</p>
</li>
<li>
<p>PRMTOP, LOADDCD<br>
<strong>Description:</strong> Paths to parameter and trajectory files.<br>
<strong>Default value:</strong> none</p>
</li>
<li>
<p>TSCALE<br>
<strong>Description:</strong> The time interval between frames in the trajectory multiplied by FSTEP.<br>
<strong>Default value:</strong> 5.0</p>
</li>
</ul>
<h3 id="optional-tfield-parameters">Optional TFIELD parameters</h3>
<ul>
<li>
<p>DENSITY<br>
<strong>Description:</strong> Grid density (1/Angstrom)<br>
<strong>Default value:</strong>  1.0</p>
</li>
<li>
<p>CUTOFF<br>
<strong>Description:</strong> Only grid cells with water occupancy higher than CUTOFF will be used for calculation of tensor field.<br>
<strong>Default value:</strong> 0.001</p>
</li>
<li>
<p>FSTEP<br>
<strong>Description:</strong> If the time interval between trajectory frames is small it is possible to increase it for computation of water displacements. For example, if FSTEP is 1 displacement is calculated from frames 1-0, 2-1, 3-2 … If FSTEP is 2, it is calculated from frames 2-0, 3-1, 4-2 …<br>
<strong>Default value:</strong>  1</p>
</li>
</ul>
<h3 id="tfield-output-options">TFIELD output options</h3>
<ul>
<li>ADC<br>
<strong>Description:</strong> Save apparent diffusion coefficient map in the occupancy field of the file ADC.pdb. The weighting factors (grid occupancy) are saved in the beta field <br>
<strong>Default value:</strong> 1 (YES)</li>
<li>FA<br>
<strong>Description:</strong> Save fractional anisotrory map in the occupancy field of the file FA.pdb<br>
<strong>Default value:</strong> 1 (YES)</li>
<li>DIFF<br>
<strong>Description:</strong> Save diffusion map in the occupancy field of the file diff.pdb<br>
<strong>Default value:</strong> 1 (YES)</li>
<li>HYDRO<br>
<strong>Description:</strong> Save water oxygen density map in the occupancy field of the file Hydrogen.pdb<br>
<strong>Default value:</strong> 1 (YES)</li>
<li>OXY<br>
<strong>Description:</strong> Save water hydrogen density map in the occupancy field of the file Oxygen.pdb<br>
<strong>Default value:</strong> 1 (YES)</li>
<li>TENSORS<br>
<strong>Description:</strong> Save diffusion tensors in the file tensors.sit<br>
<strong>Default value:</strong> 1 (YES)</li>
</ul>
<h3 id="streamline-configuration-parameters">STREAMLINE configuration parameters</h3>
<p>Run the programs:<br>
./tensors &lt; tfield.conf<br>
./streamline &lt; steamline.conf</p>
<p>Or:<br>
cd test<br>
./run_test.sh</p>
<h3 id="getting-help-with-the-configuration">Getting help with the configuration:</h3>
<p>Run “tensors” program interactively (without input from configuration files) and type “help” at the prompt.</p>
<h3 id="output-of-the-streamline-program">Output of the “streamline” program:</h3>
<ol>
<li>Streamlines color-coded by  anisotropy: streamline_A.mol2</li>
<li>Streamlines color-coded by  diffusion:  streamline_D.mol2</li>
<li>Streamlines color-coded by  direction:  streamline_XYZ.mol2</li>
</ol>
<p>The color code is saved in the charge section of mol2 files. To visualize streamlines colored by direction use shell script “load_streamlines_mol2.sh” (requires VMD).</p>
<h3 id="references">References:###</h3>
<ol>
<li>C. Gustafsson, S. Vassiliev, C. Kurten, P.-O. Syren, T. Brinck, MD simulations reveal complex water paths in squalene hopene cyclase - tunnel obstructing mutations alter the movements of water in the active site, ACS Omega. submitted (2017).</li>
<li>S. Vassiliev, P. Comte, A. Mahboob, D. Bruce, Tracking the Flow of Water through Photosystem II Using Molecular Dynamics and Streamline Tracing, Biochemistry. 49 (2010) 1873–1881. doi:10.1021/bi901900s.</li>
</ol>

