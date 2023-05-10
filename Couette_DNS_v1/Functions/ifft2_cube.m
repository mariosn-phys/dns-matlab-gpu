function [ifft2_single ] = ifft2_cube(spectral_field)
%2D Inverse Fourier Transorm
ifft2_single=real(permute(ifft2(permute(spectral_field,[2 3 1])),[3 1 2]));

end

