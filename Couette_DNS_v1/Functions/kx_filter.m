function [upert_k,upert_nk] = kx_filter(upert,k)

% filter function for streamwise harmonics of the flow field
% used to isolate streamwise harmonics, 
% k counts from 1 // k_x=1 k=2

nn=size(upert);
ny=nn(1);
NX=nn(2);
    upert_k = 0*upert;

    upert_nk = 0*upert;
            
%     for qq=1:ny
%         
%         upert_hat = fft2(squeeze(upert(qq,:,:))); 
% 
%         upert_khat = 0*upert_hat; 
% 
%         % k counts from 1 // k_x=1 k=2
%         
%         upert_khat([k NX-k+2],:) =  upert_hat([k NX-k+2],:);
%         
%         upert_k(qq,:,:) = ifft2(upert_khat); 
%         
%         upert_hat([k NX-k+2],:) = 0*upert_hat([k NX-k+2],:);
%         
%         upert_nk(qq,:,:) = ifft2(upert_hat);
%         
%     end
%     
    
    upert_hat=fft2_cube(upert);
    
    upert_khat = 0*upert_hat; 

    upert_khat(:,[k NX-k+2],:) =  upert_hat(:,[k NX-k+2],:);

    upert_k = ifft2_cube(upert_khat); 
    
    upert_hat(:,[k NX-k+2],:) = 0*upert_hat(:,[k NX-k+2],:);

    upert_nk = ifft2_cube(upert_hat);

end

