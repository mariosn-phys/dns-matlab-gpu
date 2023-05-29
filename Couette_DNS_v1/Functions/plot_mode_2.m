function plot_mode_2(kfig,un,vn,wn,gn,inx,iny,inz)
% select figure+1 and planes to plot (inx,iny,inz)
% (1) instantaneous streamwise mean (2) U velocity, plane in x
% (3) W velocity, plane in z (4) V velocity, plane in y


global yE z x

%[un,vn,wn,gn]=pick_vec(uL,vL,wL,gL,Ln,Lk);

figure(kfig+1);clf
subplot(221)

    Uyz_t=mean(un,2);
%    Uyz_t=squeeze(mean(u0,2)-repmat(UP1,[1,1,MZ]));
    Vyz_t=mean(vn,2);
    Wyz_t=mean(wn,2);

      
%    [ey,iny]=max(trapz(yE,Vyz_t.^2));

%    [ez,inz]=max(trapz(z,Wyz_t.^2,3));
   
    plot_contour(Uyz_t,Vyz_t,Wyz_t,'z','y','yz')
                title('$U(y,z,t)$','Interpreter','latex')

    subplot(223)
    
%    [ex,inx]=max(squeeze(trapz(yE,trapz(z,un.^2,3),1)));
    
    Uyz_t=un(:,inx,:);
%    Uyz_t=squeeze(mean(u0,2)-repmat(UP1,[1,1,MZ]));
    Vyz_t=vn(:,inx,:);
    Wyz_t=wn(:,inx,:);
    
    
    plot_contour(Uyz_t,Vyz_t,Wyz_t,'z','y','yz')
            title(['$u(x=',num2str(x(inx),'%0.2f'),',y,z,t)$'],'Interpreter','latex')

    subplot(222)
    
    Uyz_t=un(:,:,inz);
%    Uyz_t=squeeze(mean(u0,2)-repmat(UP1,[1,1,MZ]));
    Vyz_t=vn(:,:,inz);
    Wyz_t=wn(:,:,inz);
    
    plot_contour(Uyz_t,Vyz_t,Wyz_t,'x','y','yx')
        title(['$w(x,y,z=',num2str(z(inz),'%0.2f'),')$'],'Interpreter','latex')

    
        subplot(224)
    
    Uyz_t=un(iny,:,:);
%    Uyz_t=squeeze(mean(u0,2)-repmat(UP1,[1,1,MZ]));
    Vyz_t=vn(iny,:,:);
    Wyz_t=wn(iny,:,:);
    
    
    plot_contour(Uyz_t,Vyz_t,Wyz_t,'x','z','xz')
        title(['$v(x,y=',num2str(yE(iny),'%0.2f'),',w)$'],'Interpreter','latex')
