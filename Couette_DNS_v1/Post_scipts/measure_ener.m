
% %    Ret=sqrt(Re*Inp)
%     
% %    plot_snapshot(u0)
%     Efm(iter)=Ener(u0-Uback,v0,w0);
% %    Efp(iter,fin)=Ener(u0-repmat(mean(u0,2),[1,NX,1]),v0-repmat(mean(v0,2),[1,NX,1]),w0-repmat(mean(w0,2),[1,NX,1]));



    

%     
     O_bot(iter)=DYF(1,:)*UP1;
     O_top(iter)=DYF(end,:)*UP1;
%  
    if modf=='c'
     Inp(iter)=(O_bot(iter)+O_top(iter))/2; % Energy input Couette
    elseif modf=='p' || modf=='m'
     Inp(iter)=(O_bot(iter)-O_top(iter))/2; % Energy input Poiseuille
    end
     
     [gx,gz]=vor_xz(u0,v0,w0); % Streamwise and spanwise vorticity
 
     Dissip(iter)=2*Ener(gx,g0,gz); 
     
% Two-dimensional spectra on wall-normal planes    
% if cspec==1
%    Euu_y(:,:,:,iter)=Spec1(u0);
%    Evv_y(:,:,:,iter)=Spec1(v0);
%    Eww_y(:,:,:,iter)=Spec1(w0);
%    Egg_y(:,:,:,iter)=Spec1(g0);
%    Exx_y(:,:,:,iter)=Spec1(gx);
%    Ezz_y(:,:,:,iter)=Spec1(gz); 
%end

% Integrated Two-dimensional spectra
if cspec==1
     Euu(:,:,iter)=trapz(y,Spec1(u0),3);
     Evv(:,:,iter)=trapz(y,Spec1(v0),3);
     Eww(:,:,iter)=trapz(y,Spec1(w0),3);
     Egg(:,:,iter)=trapz(y,Spec1(g0),3);
     Exx(:,:,iter)=trapz(y,Spec1(gx),3);
     Ezz(:,:,iter)=trapz(y,Spec1(gz),3); 
end     
     
%     Time-series of the streamwise mean flow
%     vmean1=squeeze(mean(v0(:,:,:,1),2)); % vmean=0*vmean;
%     gmean1=squeeze(mean(g0(:,:,:,1),2)); % gmean=0*gmean;
%     umean1=squeeze(mean(u0(:,:,:,1),2));
%     wmean1=squeeze(mean(w0(:,:,:,1),2));
%     
%     Ust(:,:,iter)=umean1;
%     Vst(:,:,iter)=vmean1;
%     Wst(:,:,iter)=wmean1;
%     Oxt(:,:,iter)=squeeze(mean(difY_F(w0,1)-difZ_F(v0,1),2));
