# couette-dns-matlab-private
A MATLAB Script that performs Direct Numerical Simulations (DNSs) of plane Couette turbulence in a plane parallel,
streamwise and spanwise periodic 3D domain.

The provided script can be run on CPU or GPU by switching off or on the igpu flag at the beginning.

The folders contain the required functions ('Functions/'), initial fields for Re600 and Re2250 simulations ('Data/Rex/'), a script to 
acquire basic one and two point statistics ('Post_Scripts/'), 
the main executable (NL3D_Couette_gpu_on_RK3.m) 
and a script to modify the wall-normal resolution of the initial velocity field (Change_wall_normal_res.m).

It is possible to either retain the folder structure and add the base folder and subfolders to the matlab path , or merge everything in a 
single folder and edit the save and load paths accordingly.
