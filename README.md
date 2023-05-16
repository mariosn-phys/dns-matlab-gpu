# couette-dns-matlab-gpu
A MATLAB Script that performs Direct Numerical Simulations (DNSs) of plane Couette turbulence in a plane parallel,
streamwise and spanwise periodic 3D domain.

The provided script can be run on CPU or GPU by switching off or on the igpu flag at the beginning.

The folders contain the required functions ('Functions/'), initial fields for Re600 and Re2250 simulations ('Data/Rex/'), a script to 
acquire basic one and two point statistics ('Post_Scripts/'), 
the main executable (NL3D_Couette_gpu_on_RK3.m) 
and a script to modify the wall-normal resolution of the initial velocity field (Change_wall_normal_res.m).

It is possible to either retain the folder structure and add the base folder and subfolders to the matlab path , or merge everything in a 
single folder and edit the save and load paths accordingly.

Options
------------
solv | 1 to calculate and store solver matrices 0 to load stored ones (parameters must be the same!),
af | multiplies nonlinear term of perturbation-perturbation interactions, 1 for DNS,
igpu | 1 gpu is on 0 gpu is off,
mod | 'c' for Couette flow , 'p' for Poiseuille (testing purposes only), 'z' no flow,
tsav | Interval of save states,
tplot | Interval of diagnostic plots

Parameters
---------------
(These have to be the same for the stored precalculated matrices to work!)
Re | Reynolds number,
dt | time step,
a  | Fundamental Wavenumber in x,
b  | Fundamental Wavenumber in z,
N  | An odd number of grid points in y (chebyshev),
NX | An even number of grid points in x,
MZ | An even number of grid points in z
 
