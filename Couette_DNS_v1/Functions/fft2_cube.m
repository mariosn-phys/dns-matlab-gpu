function [fft2_single ] = fft2_cube(physical_field)
%2D Fourier Transform
fft2_single=permute(fft2(permute(physical_field,[2 3 1])),[3 1 2]);

end

