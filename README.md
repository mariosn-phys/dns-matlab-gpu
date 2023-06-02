# couette-dns-matlab-gpu
A MATLAB Script that performs Direct Numerical Simulations (DNSs) of plane Couette turbulence in a plane parallel,
streamwise and spanwise periodic 3D domain.

The provided script can be run on CPU or GPU by switching off or on the igpu flag at the beginning.

The folders contain the required functions ('Functions/'), initial fields for Re600, Re2250 and R3000(coarse grid) simulations ('Data/Rex/'), a script to 
acquire basic one and two point statistics ('Post_Scripts/Avg_statistics_plot.m'), 
the main executable (Scripts/NL3D_Couette_gpu_on_RK3.m) 
and a script to modify the wall-normal resolution of the initial velocity field (Scripts/Change_wall_normal_res.m). A folder 'Files/' should be created to output
the modified initial fields (or select a different path).

It is possible to either retain the folder structure and add the base folder and subfolders to the matlab path , or merge everything in a 
single folder and edit the save and load paths accordingly.

Options
------------
field_path | set path of DNS save state folder,
fmt | format of time in filenames, default is '%04.2f'
diag_file | filename for basic diagnostic quantities (Input, energy, CFL)
Ti | Starting time  (must match an existing state in the field_path),
Tf | Time at the last integration step,
solv | 1 to calculate and store solver matrices 0 to load stored ones (parameters must be the same!),
psolv | 1 enables parallel pool to build solvers, 0 for serial mode
npc | number of cores to utilize in the parallel pool
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
N  | An odd number of grid points in y (chebyshev),
NX | An even number of grid points in x,
MZ | An even number of grid points in z
 
(Sets of parameters can be found in the initial field folders 'Data/Rex')
