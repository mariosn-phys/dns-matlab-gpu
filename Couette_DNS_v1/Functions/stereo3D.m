function stereo3D(u0,v0,w0,fig,str)

global a b


% Resolution

NN=size(u0);
%cut=0.1;
sh=0;


% Grid
N=NN(1)-2;NX=NN(2);MZ=NN(3);


L=2;dy=L/(N+1);y=linspace(-1+dy,1-dy,N)';yE=[-1;y;1];
A=2*pi/a;dx=A/NX;x=-dx+dx*(1:NX)';xE=[x;A]; 
B=2*pi/b;dz=B/MZ;z=-dz+dz*(1:MZ)';zE=[z;B];

%[ym,xm,zm]=meshgrid(x,yE,z);
[zm,xm,ym]=ndgrid(z,x,yE);

if str == 0
u0=u0-repmat(mean(u0,2),[1 NX 1]);
v0=v0-repmat(mean(v0,2),[1 NX 1]);
w0=w0-repmat(mean(w0,2),[1 NX 1]);   
else    
u0=u0-repmat(mean(mean(u0,3),2),[1 NX MZ]);
v0=v0-repmat(mean(mean(v0,3),2),[1 NX MZ]);
w0=w0-repmat(mean(mean(w0,3),2),[1 NX MZ]);
end

% u2=squeeze(mean(u0.^2,2));
% v2=squeeze(mean(v0.^2,2));
% w2=squeeze(mean(w0.^2,2));
% 
% figure(3+fig);clf;
% contourf(z,yE,u2)
% 
% figure(4+fig);clf;
% contourf(z,yE,v2)
% 
% figure(5+fig);clf;
% contourf(z,yE,w2)

u0=circshift(u0,[0 0 sh]);
u0=permute(repmat(mean(u0,2),[1 NX 1]),[3 2 1]);


v0=circshift(v0,[0 0 sh]);
v0=permute(v0,[3 2 1]);


w0=circshift(w0,[0 0 sh]);
w0=permute(w0,[3 2 1]);

cf=max(u0(:));cut=0.5*cf;


figure(3+fig);clf;
%subplot(211)
for ij=1:1
isosurface(xm,zm,ym,u0,ij*cut)
hold on
isosurface(xm,zm,ym,u0,-ij*cut)
end
title(['u_max=',num2str(max(u0(:)))])
ylim([-1 1])
axis equal
    ylabel('$z$','Interpreter','Latex','Fontsize',20)
    zlabel('$y$','Interpreter','Latex','Fontsize',20)
    xlabel('$x$','Interpreter','Latex','Fontsize',20)
    set(gca,'Fontsize',20)

figure(4+fig);clf
%subplot(212)
for ij=1:1
isosurface(xm,zm,ym,v0,ij*0.5*cut)
hold on
isosurface(xm,zm,ym,v0,-ij*0.5*cut)
end
title(['v_max=',num2str(max(u0(:)))])
axis equal

figure(5+fig);clf
% contourf
for ij=1:1
isosurface(xm,zm,ym,w0,ij*0.5*cut)
hold on
isosurface(xm,zm,ym,w0,-ij*0.5*cut)
end
title(['w_max=',num2str(max(u0(:)))])
axis equal

end
%isosurface(u0,-0.1)
%hold on
%isosurface(u0,0.1)