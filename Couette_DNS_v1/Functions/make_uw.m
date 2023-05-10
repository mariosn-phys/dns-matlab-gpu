function [u0,v0,w0,g0] = make_uw(vi,gi,UP1,WP1,b1,b2)
%Employing incompressibility condition to reconstruct the physical
%flow field from v and g

global NX MZ k l DYF

for kk=1:NX/3+1
    
    kp=k(kk);
    
    for jj=[1:MZ/3+1 MZ-MZ/3+1:MZ]
        
        lp=l(jj);
        
        vihat=vi(:,kk,jj);
        gihat=gi(:,kk,jj);
        
        if (kp==0) && (lp==0)
        else
            dydvihat=DYF*vihat;
            ui(:,kk,jj)=(-1i*lp*gihat+1i*kp*dydvihat)/(kp^2+lp^2);
            wi(:,kk,jj)=(1i*kp*gihat+1i*lp*dydvihat)/(kp^2+lp^2);
        end
        
    end

end

for kk=2:NX/3+1
    
    gi(:,NX-kk+2,1)=conj(gi(:,kk,1));
    vi(:,NX-kk+2,1)=conj(vi(:,kk,1));
    ui(:,NX-kk+2,1)=conj(ui(:,kk,1));
    wi(:,NX-kk+2,1)=conj(wi(:,kk,1));
    
    jj=[2:MZ/3+1 MZ-MZ/3+1:MZ];
    gi(:,NX-kk+2,MZ-jj+2)=conj(gi(:,kk,jj));
    vi(:,NX-kk+2,MZ-jj+2)=conj(vi(:,kk,jj));
    ui(:,NX-kk+2,MZ-jj+2)=conj(ui(:,kk,jj));
    wi(:,NX-kk+2,MZ-jj+2)=conj(wi(:,kk,jj));
        
end

v0=ifft2_cube(vi);
u0=ifft2_cube(ui);
w0=ifft2_cube(wi);
g0=ifft2_cube(gi);

u0=repmat(UP1,[1,NX,MZ])+u0;u0=[ones(1,NX,MZ)*b1;u0(2:end-1,:,:);ones(1,NX,MZ)*b2];
w0=repmat(WP1,[1,NX,MZ])+w0;w0=[zeros(1,NX,MZ);w0(2:end-1,:,:);zeros(1,NX,MZ)];
