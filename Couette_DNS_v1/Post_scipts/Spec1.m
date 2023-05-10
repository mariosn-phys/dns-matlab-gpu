function Euu_total = Spec1(U)

%%%% Adapted from the Channel Code version written by Adrian L.Duran %%%% 


global N NX MZ dy dx dz A B L 

ny=N;

U_hat = fft2(permute(U(2:end-1,:,:),[2 3 1])); %changed ifft to fft

% dealias
mx    = 2/3*NX;
i1    = NX/2+1-NX/3/2;
i2    = NX/2+1+NX/3/2+1;
U_hat = U_hat([1:i1 i2:NX],:,:);

mz    = 2/3*MZ;
k1    = MZ/2+1-MZ/3/2;
k2    = MZ/2+1+MZ/3/2+1;
U_hat = U_hat(:,[1:k1 k2:MZ],:);

% getting spectra
ii  = [1:mx/2 mx/2 mx/2:-1:2];
kk  = [1:mz/2 mz/2 mz/2:-1:2];

Euu = zeros(mx/2,mz/2,ny);
for i=1:mx
   for k=1:mz
       Euu(ii(i),kk(k),:) = Euu(ii(i),kk(k),:) + abs(U_hat(i,k,:)).^2; %SPA1
   end
end

Euu=Euu/NX/MZ;

Euu_total=Euu/2*dx*dz/A/B/L;

end
