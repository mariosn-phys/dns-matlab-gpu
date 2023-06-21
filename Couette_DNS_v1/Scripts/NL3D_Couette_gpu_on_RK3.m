
%    A DNS Script of plane parallel flows (Couette/Poiseuille) for MATLAB   
%    Copyright (C) 2023 Marios-Andreas Nikolaidis

%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU Affero General Public License as published
%    by the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU Affero General Public License for more details.

%    You should have received a copy of the GNU Affero General Public License
%    along with this program.  If not, see <https://www.gnu.org/licenses/>.

% Base DNS code used in Nikolaidis & Ioannou (2022)
% A pseudo-spectal Navier-Stokes solver for plane parallel Couette flow
% in wall-normal velocity/ vorticity formulation (Kim,Moin,Moser 1987 - 
% an example of a spectral mpi DNS has been developed by the Fluid Dynamics 
% group at UPM/ J.Jimenez et al.), that is capable of simulating turbulent flow
% dynamics.
% Differentiations are performed with the pseudo-spectral Chebyshev
% (Weideman & Reddy 2000) and Fourier (e.g Trefethen 2000) matrices. 
% Tested with NVIDIA gpus, enabled with igpu flag. Exceeding available GPU memory
% will result in crashes.
% Wherever possible during the time-stepping for-loops have been restructured to 
% matrix multiplication operations, which utilizes the build-in vectorization of 
% MATLAB (see also the DNS codes of Vuorinen & Keskinen 2015)

	
	
	
	
clear all;

global M N Re NX D1x D2x D1z D2z DYF D2F D2 D4 y gamma
global MZ k l dt xE yE zE A B L dx dz kkm llm DYkron 
 
global ICvkron1 ICvDvkron1 ICgkron1 ICggkron1
global ICvkron2 ICvDvkron2 ICgkron2 ICggkron2

% global ICvkron1c ICvDvkron1c ICgkron1c ICggkron1c
% global ICvkron2c ICvDvkron2c ICgkron2c ICggkron2c
% 
% global ICvCvRK1 ICvIDELvRK1 ICgRK1 ICgCgRK1
% global ICvCvRK2 ICvIDELvRK2 ICgRK2 ICgCgRK2

global S_mf S_mp Sol_m b1 b2
% global kLyap NLyap kkmm llmm DYkronm
global a b

%Function paths or add main folder ('../') and subfolders to path
addpath('../Functions/')
%addpath('Functions_modal/')

solv=1; % 1 to calculate and save solver matrices, 0 to load stored 
%precalculated ones for the same parameters
psolv=0; % 1 enables parallel pool for solver matrices, 0 serial mode
npc=2; % number of parallel workers in pool for solvers, increases 
       %memory requirements during calculation

af=1;  
% Nonlinear term modulation, 0 RNL 1 DNS // If used only as DNS 
% advect_NL_anl is changed to advect_NL for faster execution

%%% GPU ---------------------------------
igpu=0; % 1 gpu on, 0 gpu off
if igpu
gpuDevice(1) % Assign to gpudevice
end
%%% -------------------------------------

%%% !!!!!!!!!!!!
field_path = '../Data/Re3000_n71_coarse/'  %% save and restart file path
diag_file  = 'diagnostics'; %% save elementary diagnostics
%%% !!!!!!!!!!!!

Re=3000; % Reynolds number 

% Poiseuille (test) 'p' or Couette 'c'
modf='c';
it0=2; % First advance step
dt=0.008; % Plot and save intervals assume that 1/dt = integer. Some compatible 
         % choices for dt=[0.025,0.02,0.0125,0.01,0.008,0.00625,0.005] 
Ti=5499; % Initial Time (should match time on restart file name)
Tf=5700; % Final Time
T=Ti:dt:Tf; 
NT=length(T);
g=dt/(2*Re); % CN coefficient

%%% File name %% Save file defined at the end of the time-stepping loop
fmt='%04.2f'; % Load and save time label format
start_file=[field_path,'state_Re',num2str(Re),'_',num2str(T(1),fmt),'.mat']

tsav=1; % interval of saves
tplot=100; % plot basic diagnostics

% Fundamental wavenumbers
%Re3000 and Re3250
a=2/2;
b=2/1;

%Re600 and Re2250
%a=2/1.75;
%b=2/(1.2);


%%%% NX,MZ 12,24,48 preferably increment at steps of 12 point, also works
%%%% with increments of 6. N is odd for centerline point.
%%%% Boundary points not counted. 
	
NX=108;
MZ=108;
N=71;

A=2*pi/a;dx=A/NX;x=-dx+dx*(1:NX)';xE=[x;A]; 
B=2*pi/b;dz=B/MZ;z=-dz+dz*(1:MZ)';zE=[z;B];
L=2;


%save([field_path,'parameters_Re',num2str(Re),'_n',num2str(N),'.mat'],'a','b','Re','N','NX','MZ','dt')
%Parameters for 600, 2250 and 3000 can be found in their respective folders '../Data/Re*/

% Form Differentiation matrices
% x matrix (see eg Trefethen 2000)
    column=[0 .5*(-1).^(1:NX-1).*cot((1:NX-1)*a*dx/2)];
    D1x=a*toeplitz(column,column([1 NX:-1:2]));
    D2x=D1x^2;
% z matrix    
    column=[0 .5*(-1).^(1:MZ-1).*cot((1:MZ-1)*b*dz/2)];
    D1z=b*toeplitz(column,column([1 MZ:-1:2]));
    D2z=D1z^2;
  
% y matrix (see Weideman and Reddy 2000)
    DYF=cheb(N+1);DYF=flip(flip(DYF,1),2);
    DY=DYF(2:end-1,2:end-1);
    DYkron=kron(speye((NX/3)*(2*MZ/3-1)),DY);

    D2F=DYF^2;
    D2=D2F(2:end-1,2:end-1);

    [y,D4]=cheb4c(N+2);D4=flip(flip(D4,1),2);
    y=-y;
    yE=[-1;y;1];
    
%-------%
%Required only to plot the instantaneous flow fields
%[xx,yy]=meshgrid(x,y);
%[xxE,yyE]=meshgrid(x,yE);
%-------%
    
%

M=[0:(NX/2-1) 0 (1-NX/2):(-1)];
    k=2*pi*M/A;
P=[0:(MZ/2-1) 0 (1-MZ/2):(-1)];
    l=2*pi*P/B;  
  
    if modf=='p'
    U1=4/3*(1-yE.^2);gm=-8/3/Re; %Poiseuille (constant pressure) 
    elseif modf=='c'
    U1=yE;gm=0; %Couette
    elseif modf=='z'
    U1=yE*0;gm=-2/Re*0;    
    end
    b1=U1(1);b2=U1(end);
  
    Uback=repmat(U1,[1,NX,MZ]);
    gamma=gm*ones(N+2,1);

%----------Build Solvers  
if solv==1
    if psolv==1    
    parpool(npc) 

    [ICvkron1,ICvDvkron1,ICgkron1,ICggkron1,~,~,~,~,S_mf(:,:,1),S_mp(:,:,1),Sol_m(:,:,1),kkm,llm]=parsolvers(g,1/2);
      
    [ICvkron2,ICvDvkron2,ICgkron2,ICggkron2,~,~,~,~,S_mf(:,:,2),S_mp(:,:,2),Sol_m(:,:,2),~,~]=parsolvers(g,1);

    delete(gcp)
    else
  
    [ICvkron1,ICvDvkron1,ICgkron1,ICggkron1,~,~,~,~,S_mf(:,:,1),S_mp(:,:,1),Sol_m(:,:,1),kkm,llm]=solvers(g,1/2);
      
    [ICvkron2,ICvDvkron2,ICgkron2,ICggkron2,~,~,~,~,S_mf(:,:,2),S_mp(:,:,2),Sol_m(:,:,2),~,~]=solvers(g,1);   
    end

save([field_path,'Solvers.mat'],'ICvkron1','ICvkron2','ICvDvkron1','ICvDvkron2','ICgkron1','ICgkron2','ICggkron1','ICggkron2','S_mf','S_mp','Sol_m','kkm','llm','-v7.3')
else 
%----------Load Solvers (only if calculated with same parameteres!!!)
load([field_path,'Solvers.mat'])
display('Solvers loaded')
end

display('Solver initialization complete')

% Transfer Variables to GPU 
if igpu==1
      ICvkron1=gpuArray(ICvkron1);
      ICvDvkron1=gpuArray(ICvDvkron1);
      ICgkron1=gpuArray(ICgkron1);
      ICggkron1=gpuArray(ICggkron1);
      kkm=gpuArray(kkm);
      llm=gpuArray(llm);
      DYkron=gpuArray(DYkron);

      ICvkron2=gpuArray(ICvkron2);
      ICvDvkron2=gpuArray(ICvDvkron2);
      ICgkron2=gpuArray(ICgkron2);
      ICggkron2=gpuArray(ICggkron2);
      display('GPU upload complete')
end
%----------Build Solvers End   

%    Load init and transform to physical space
 [vi,gi,UP1,WP1]=read_from_disk_compact(start_file);
 [u0,v0,w0,g0] = make_uw(gather(vi),gather(gi),gather(UP1),gather(WP1),b1,b2);
  
% Mean flow        
%    vmean1=repmat(mean(v0(:,:,:,1),2),[1 NX 1 1]); % vmean=0*vmean;
%    gmean1=repmat(mean(g0(:,:,:,1),2),[1 NX 1 1]); % gmean=0*gmean;
%    umean1=repmat(mean(u0(:,:,:,1),2),[1 NX 1 1]);
%    wmean1=repmat(mean(w0(:,:,:,1),2),[1 NX 1 1]);
    
%     figure(91);
%     pcolor(squeeze(umean1(:,1,:)))
    
    if igpu==1
         v0=gpuArray(v0);
         u0=gpuArray(u0);
         w0=gpuArray(w0);
         g0=gpuArray(g0);
    end
    
    Efm=zeros(NT,1);CFL=Efm;O_bot=Efm;O_top=Efm;
             
    Efm(1)=gather(Ener(u0-Uback,v0,w0)); % Energy of deviations from the laminar flow
    CFL(1)=gather((max(u0(:))/dx+max(abs(max(v0(2:end,:),[],2)./diff(yE)))+max(w0(:))/dz)*dt);
	
    % Shear at top and bottom walls	
    O_bot(1)=gather(DYF(1,:)*mean(mean(u0,3),2)); 
    O_top(1)=gather(DYF(end,:)*mean(mean(u0,3),2));
 
    
    %%
    tic;
    tdt=[];
for it=it0:NT
    
   %if rem(it,125)==1
   tic;
   %end
   
    % u init
    
    up=u0;
    vp=v0;
    wp=w0;
    gp=g0;
    
    if af==1
    [advv1,advg1] = advect_NL(up,vp,wp,gp);   
    else
    [advv1,advg1] = advect_NL_anl(up,vp,wp,gp,af);
    end
    
    advU1=DYF*mean(mean(up.*vp,3),2)+gamma;
    advW1=DYF*mean(mean(wp.*vp,3),2);
    
    [ui,vi,wi,gi,UP1,WP1] = rkstep(u0,v0,w0,g0,advg1,advv1,advU1,advW1,1,1/2);
    
    v1=ifft2_cube(vi);
    u1=ifft2_cube(ui);
    w1=ifft2_cube(wi);
    g1=ifft2_cube(gi);

    u1=repmat(UP1,[1,NX,MZ])+u1;u1=[ones(1,NX,MZ)*b1;u1(2:end-1,:,:);ones(1,NX,MZ)*b2];
    w1=repmat(WP1,[1,NX,MZ])+w1;w1=[zeros(1,NX,MZ);w1(2:end-1,:,:);zeros(1,NX,MZ)];
    
    % u step 1
    
    up=u1;
    vp=v1;
    wp=w1;
    gp=g1;
    
    if af==1
    [advv2,advg2] = advect_NL(up,vp,wp,gp);   
    else
    [advv2,advg2] = advect_NL_anl(up,vp,wp,gp,af);
    end 
        
    advU2=DYF*mean(mean(up.*vp,3),2)+gamma;
    advW2=DYF*mean(mean(wp.*vp,3),2);
        
    [ui,vi,wi,gi,UP2,WP2] = rkstep(u0,v0,w0,g0,2*advg2-advg1,2*advv2-advv1,2*advU2-advU1,2*advW2-advW1,2,1);

    v2=ifft2_cube(vi);
    u2=ifft2_cube(ui);
    w2=ifft2_cube(wi);
    g2=ifft2_cube(gi);

    u2=repmat(UP2,[1,NX,MZ])+u2;u2=[ones(1,NX,MZ)*b1;u2(2:end-1,:,:);ones(1,NX,MZ)*b2];
    w2=repmat(WP2,[1,NX,MZ])+w2;w2=[zeros(1,NX,MZ);w2(2:end-1,:,:);zeros(1,NX,MZ)];
    
    % u step 2
    
    up=u2;
    vp=v2;
    wp=w2;
    gp=g2;
    
    if af==1
    [advv3,advg3] = advect_NL(up,vp,wp,gp);   
    else
    [advv3,advg3] = advect_NL_anl(up,vp,wp,gp,af);
    end
    
    advU3=DYF*mean(mean(up.*vp,3),2)+gamma;
    advW3=DYF*mean(mean(wp.*vp,3),2);
        
    [ui,vi,wi,gi,UP3,WP3] = rkstep(u0,v0,w0,g0,(advg1+4*advg2+advg3)/6,(advv1+4*advv2+advv3)/6,(advU1+4*advU2+advU3)/6,(advW1+4*advW2+advW3)/6,2,1);

    v0=ifft2_cube(vi);
    u0=ifft2_cube(ui);
    w0=ifft2_cube(wi);
    g0=ifft2_cube(gi);

    % u final
    
    u0=repmat(UP3,[1,NX,MZ])+u0;u0=[ones(1,NX,MZ)*b1;u0(2:end-1,:,:);ones(1,NX,MZ)*b2];
    w0=repmat(WP3,[1,NX,MZ])+w0;w0=[zeros(1,NX,MZ);w0(2:end-1,:,:);zeros(1,NX,MZ)];  
            
%     [~,u0]=kx_filter(u0,srg);
%     [~,v0]=kx_filter(v0,srg);
%     [~,w0]=kx_filter(w0,srg);
%     [~,g0]=kx_filter(g0,srg);

    Efm(it)=gather(Ener(u0-Uback,v0,w0)); 
    
    CFL(it)=gather((max(u0(:))/dx+max(abs(max(v0(2:end,:),[],2)./diff(yE)))+max(w0(:))/dz)*dt);
   
    O_bot(it)=gather(DYF(1,:)*mean(mean(u0,3),2));
    O_top(it)=gather(DYF(end,:)*mean(mean(u0,3),2));

    if rem(T(it),tplot)==0
    figure(99);
    plot(UP3,yE,'+',U1,yE,'r')
    xlabel('U'),ylabel('y')
    figure(105);
    contourf(x,yE,v0(:,:,1),20)
    xlabel('x');ylabel('y');title('v')
    figure(111);clf
    plot(T(1:it),Efm(1:it))
    xlabel('T');ylabel('E')
    figure(112);clf
    plot(T(1:it),CFL(1:it))
    xlabel('T');ylabel('CFL')
    figure(113);clf
    	if modf='c'
    	plot(T(1:it),sqrt(Re*(O_bot(1:it)+O_top(1:it))/2))
    	elseif modf='p'
    	plot(T(1:it),sqrt(Re*(O_bot(1:it)-O_top(1:it))/2))
	end
    xlabel('T');ylabel('Re_{\tau}')
    figure(122);
    contourf(z,yE,squeeze(mean(u0,2)))
    xlabel('z');ylabel('y');title('U(y,z,t)')
    drawnow
    figure(123);
    contourf(z,yE,squeeze(mean(u0-repmat(mean(mean(u0,3),2),[1 NX MZ]),2)))
    xlabel('z');ylabel('y');title('streak')
    drawnow
    end
    

   
    %if rem(it,125)==0
    tdt=[tdt toc]; % time step bench
    %end
     
    if rem(T(it),tsav)==0
    
    % Write Fourier coefficients to disk    
    write_to_disk_compact(gather(vi),gather(gi),gather(UP3),gather(WP3),T(it),[field_path,'state_Re',num2str(Re),'_',num2str(T(it),fmt),'.mat']);
    
    save([field_path,diag_file,'.mat'],'Efm','CFL','O_bot','O_top','tdt','T')

    end
    
    % Code Breakdown check
    %if max(abs(u0(:)))>=1.2*abs(max(U1))
    %   display('error')
    %	break
    %end
    
end
    
