function plot_contour(u2,v2,w2,xlb,ylb,ori)
% contour cross-sections of the velocity normal
% to the selected ori - plane

global x yE z
global xx_yx yy_yx yy_yz zz_yz xx_xz zz_xz

umax=max(abs(u2(:)));u2=u2/umax;

if ori=='yz'

    Uyz_t=squeeze(u2);
%    Uyz_t=squeeze(mean(u0,2)-repmat(UP1,[1,1,MZ]));
    Vyz_t=squeeze(v2);
    Wyz_t=squeeze(w2);
       
 %   figure(2);clf;
    hold on
    contourf(z,yE,Uyz_t);colorbar;
    quiver(zz_yz,yy_yz,Wyz_t',Vyz_t',1.5)
    xlabel(xlb);ylabel(ylb);
    xlim([0 z(end)])
    hold off

elseif ori=='yx'
    
    Uyz_t=squeeze(u2);
%    Uyz_t=squeeze(mean(u0,2)-repmat(UP1,[1,1,MZ]));
    Vyz_t=squeeze(v2);
    Wyz_t=squeeze(w2);
       
 %   figure(2);clf;
    hold on
    contourf(x,yE,Wyz_t);colorbar;
    quiver(xx_yx,yy_yx,Uyz_t,Vyz_t,1.5)
    xlabel(xlb);ylabel(ylb);
    xlim([0 x(end)])
    hold off
    
elseif ori=='xz'
    
    Uyz_t=squeeze(u2);
%    Uyz_t=squeeze(mean(u0,2)-repmat(UP1,[1,1,MZ]));
    Vyz_t=squeeze(v2);
    Wyz_t=squeeze(w2);
       
 %   figure(2);clf;
    hold on
    contourf(x,z,Vyz_t');colorbar;
    quiver(xx_xz,zz_xz,Uyz_t,Wyz_t,1.5)
    xlabel(xlb);ylabel(ylb);
    xlim([0 x(end)]);ylim([0 z(end)])
    hold off  
    
end
colormap(jet)
set(gca,'LineWidth',2,'FontSize',14)
end

