
    O_bot(iter)=DYF(1,:)*mean(mean(u0,3),2);
    O_top(iter)=DYF(end,:)*mean(mean(u0,3),2);

    Inp(iter)=(O_bot(iter)+O_top(iter))/2;

    [gx,gz]=vor_xz(u0,v0,w0);

    Dissip(iter)=2*Ener(gx,g0,gz);
%    Ret=sqrt(Re/2*(O_bot+O_top))
    
%    plot_snapshot(u0)
    Efm(iter)=Ener(u0-Uback,v0,w0);
%    Efp(iter,fin)=Ener(u0-repmat(mean(u0,2),[1,NX,1]),v0-repmat(mean(v0,2),[1,NX,1]),w0-repmat(mean(w0,2),[1,NX,1]));

%     Euu(:,:,:,iter)=Spec1(u0);
%     Evv(:,:,:,iter)=Spec1(v0);
%     Eww(:,:,:,iter)=Spec1(w0);
%     Egg(:,:,:,iter)=Spec1(g0);
%     Exx(:,:,:,iter)=Spec1(gx);
%     Ezz(:,:,:,iter)=Spec1(gz); 

    Euu(:,:,iter)=trapz(y,Spec1(u0),3);
    Evv(:,:,iter)=trapz(y,Spec1(v0),3);
    Eww(:,:,iter)=trapz(y,Spec1(w0),3);
    Egg(:,:,iter)=trapz(y,Spec1(g0),3);
    Exx(:,:,iter)=trapz(y,Spec1(gx),3);
    Ezz(:,:,iter)=trapz(y,Spec1(gz),3); 

    


