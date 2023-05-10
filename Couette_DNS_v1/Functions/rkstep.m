function [ui,vi,wi,gi,UP1,WP1] = rkstep(u0,v0,w0,g0,advg1,advv1,advU1,advW1,step,h)

global S_mf S_mp Sol_m b1 b2 DYF dt gamma


% u0 the base flow
% u1 the advection flow
          
    UP=mean(mean(u0,3),2);
    WP=mean(mean(w0,3),2);

    [gi,vi,ui,wi] = solv_vg_f_kron_zx_RK3(g0,v0,advg1,advv1,step,h);
    
    dUP=S_mp(:,:,step)*UP-h*dt*advU1;
    dWP=S_mp(:,:,step)*WP-h*dt*advW1;
    
    UP1=[b1;Sol_m(:,:,step)*(dUP(2:end-1)-S_mf(2:end-1,1,step)*b1-S_mf(2:end-1,end,step)*b2);b2];  
    WP1=[0;Sol_m(:,:,step)*dWP(2:end-1);0];
        
end

