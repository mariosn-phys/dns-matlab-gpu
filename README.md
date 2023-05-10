# couette-dns-matlab-private
A MATLAB Script that performs Direct Numerical Simulations (DNSs) of plane Couette turbulence in 3D domains. 
The provided script can be run on CPU or GPU by switching off or on the igpu flag at the beginning.

The folders contain the required functions, initial fields for Re600 and Re2250 simulations, a script to 
acquire basic one and two point statistics, the main executable and a script to modify the wall-normal 
resolution of the initial velocity field.

It is possible to either add the base folder and subfolders to the matlab path, or merge everything in a 
single folder.

%Functions folder should be added to the path before running the script
