# dns-matlab-gpu
A MATLAB Script that performs Direct Numerical Simulations (DNSs) of plane Couette turbulence in a plane parallel,
streamwise and spanwise periodic 3D domain.

The main executable files are in the folder 'Scripts/':
The NL3D_Couette_gpu_on_RK3.m script version supports simulations of Couette and Poiseuille ( constant pressure gradient)
Set up with a coarse grid Couette flow case at R=3000.

The NL3D_Poiseuille_gpu_on_RK3.m script version supports simulations of Couette and Poiseuille ( constant pressure gradient or mass flux)
Set up with a coarse grid Poiseuille flow case at R=3250.

The provided scripts run on CPU or GPU by switching off or on the igpu flag at the beginning.

Initial fields for different cases are provided in the folders found in 'Data/'. The parameters of these simulations can be found in the 'parameters.mat' file
in each folder. 

Different grid resolution in x and z is automatically adjusted when loading a file. To modify the wall-normal resolution of the initial velocity field, 
an interpolation script is provided (Scripts/Change_wall_normal_res.m). A folder 'Files/' should be created to output the modified initial fields 
(or select a different path).

Basic one- and two- point statistics can be acquired with the 'Avg_statistics_plot.m' script found in the 'Post_Scripts/' folder. This script calls 
a function for plotting slices of the velocity field. It is possible to calculates the spectra of wall-normal planes by commenting the appropriate 
section in the 'measure_plot_ener.m' script.

The required functions to run these scripts are found in the folder ('Functions/'). It is possible to either retain the folder structure and add the base folder 
and subfolders to the matlab path , or merge everything in a single folder and edit the save and load paths accordingly.

Options
------------
field_path | set path of DNS save state folder,
fmt | format of time in filenames, default is '%04.2f',
diag_file | filename for basic diagnostic quantities (Input, energy, CFL),
Ti | Starting time  (must match an existing state in the field_path),
Tf | Time at the last integration step,
solv | 1 to calculate and store solver matrices 0 to load stored ones (parameters must be the same!),
psolv | 1 enables parallel pool to build solvers, 0 for serial mode,
npc | number of cores to utilize in the parallel pool,
af | multiplies nonlinear term of perturbation-perturbation interactions in the perturbation equation, 1 for DNS,
igpu | 1 gpu is on 0 gpu is off,
modf | 'c' for Couette flow , 'p' for Poiseuille constant pressure, 'm' for Poiseuille constant mass, 'z' no flow,
tsav | Interval of save states (product of integer times dt),
tplot | Interval of diagnostic plots (product of integer times dt)

Parameters
---------------
(When solver matrices are loaded these parameters have to be the same with the ones used to precalculate them in order for the DNS to work correctly!)
Re | Reynolds number,
dt | time step,
a  | Fundamental Wavenumber in x,
b  | Fundamental Wavenumber in z,
N  | An odd number of grid points in y (Chebyshev),
NX | An even number of grid points in x,
MZ | An even number of grid points in z
 
(Sets of parameters can be found in the initial field folders 'Data/Rex')
