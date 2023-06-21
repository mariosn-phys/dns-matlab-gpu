clear all;

global M N Re NX D1x D2x D1z D2z DYF D2F D2 D4 y gamma
global MZ k l dt xE yE zE A B L dx dz
global b1 b2
global a b kkm llm kkmm llmm DYkron DYkronm
global x z
global xx_yx yy_yx yy_yz zz_yz xx_xz zz_xz

%gpuDevice(1)


addpath('../Functions/')

%%%% Options
field_path='../Data/Re3250_n81/'  %% save and restart file path
fname=['P3250T']; %% Save file name
fmt='%04.2f'; %% Format of time on filenames
cspec=0; % Calculate 2d spectra 
% Poiseuille 'm' 'p' or Couette 'c'
modf='m';

Re=3250;%Re=dum(2);
Ti=4600;Tf=4630;Tstp=1;
T=Ti:Tstp:Tf;
NT=length(T);

NX=144;
MZ=144;
N=81;

%load([field_path,'parameters_Re',num2str(Re),'_n',num2str(N),'.mat'],'a','b','Re','N','NX','MZ')

% a=2/(1.75);%a=2*pi/1.774;
% b=2/(1.2);%b=2*pi/0.884;
% 
a=2/2;
b=2/1;

%%%% NX,MZ 12,24,48, +12 or +6

A=2*pi/a;dx=A/NX;x=-dx+dx*(1:NX)';xE=[x;A]; 
B=2*pi/b;dz=B/MZ;z=-dz+dz*(1:MZ)';zE=[z;B];
L=2;




% x matrix (Trefethen 2000)
    column=[0 .5*(-1).^(1:NX-1).*cot((1:NX-1)*a*dx/2)];
    D1x=a*toeplitz(column,column([1 NX:-1:2]));
    D2x=D1x^2;
% z matrix    
    column=[0 .5*(-1).^(1:MZ-1).*cot((1:MZ-1)*b*dz/2)];
    D1z=b*toeplitz(column,column([1 MZ:-1:2]));
    D2z=D1z^2;
  
% y matrix (Weideman and Reddy 2000)
    DYF=cheb(N+1);DYF=flip(flip(DYF,1),2);
    DY=DYF(2:end-1,2:end-1);
    DYkron=kron(speye((NX/3)*(2*MZ/3-1)),DY);
    DYkronm=kron(speye(sum(NLyap/2)*(2*MZ/3-1)),DY);
    
    
    D2F=DYF^2;
    D2=D2F(2:end-1,2:end-1);

    [y,D4]=cheb4c(N+2);D4=flip(flip(D4,1),2);
    y=-y;
    yE=[-1;y;1];
    
%-------% Mesh grids for plots
[xx_yx,yy_yx]=meshgrid(x,yE);
[yy_yz,zz_yz]=meshgrid(yE,z);
[xx_xz,zz_xz]=meshgrid(x,z);
%-------%

M=[0:(NX/2-1) 0 (1-NX/2):(-1)];
    k=2*pi*M/A; 
P=[0:(MZ/2-1) 0 (1-MZ/2):(-1)];
    l=2*pi*P/B;  
  
    if modf=='c'
    U1=yE;gm=0;pm=0;pp=0; display('Couette');
    elseif modf=='m'
    U1=4/3*(1-yE.^2);gm=1/2/Re*(DYF(end,:)*U1-DYF(1,:)*U1);pm=1;pp=0; display('Poiseuille (constant mass flux)'); 
    elseif modf=='p'
    U1=4/3*(1-yE.^2);gm=1/2/Re*(DYF(end,:)*U1-DYF(1,:)*U1);pm=0;pp=1; display('Poiseuille (constant pressure)');
    elseif modf=='z'
    U1=yE*0;gm=-2/Re*0;pm=0;pp=0;    
    end
    b1=U1(1);b2=U1(end);

    %Initial Conditions vortex+optimal n_x=1,n_z=1
  
    Uback=repmat(U1,[1,NX,MZ]);
    gamma=gm*ones(N+2,1);

  
   iter=0;
   
   yp=21; % y cross section to plot
%%
    for it=1:NT
    iter=iter+1
    
    start_file=[field_path,'state_Re',num2str(Re),'_',num2str(T(it),fmt),'.mat']

    [vi,gi,UP1,WP1]=read_from_disk_compact(start_file);
    [u0,v0,w0,g0] = make_uw(gather(vi),gather(gi),gather(UP1),gather(WP1),b1,b2);

    upvp(:,iter)=mean(mean(u0.*v0,3),2);
    vpvp(:,iter)=mean(mean(v0.*v0,3),2);
    wpwp(:,iter)=mean(mean(w0.*w0,3),2);
    upup(:,iter)=mean(mean((u0-repmat(gather(UP1),[1 NX MZ])).^2,3),2);

    UP(:,iter)=gather(UP1);
   
    if rem(T(it),100)==0
    plot_mode_2(17,u0,u0,u0,g0,22,32,68,'no') % Plot velocity state
    drawnow
    end
    
    measure_ener %% Input, Dissipation and Spectra

       
    end
    
 
    rg=1:iter; % range of fields used to calculate statistics
    Ret=sqrt(Inp*Re); % Friction Reynolds number


save([field_path,'statistics_uu_',fname,'.mat'],'upvp','UP','upup','vpvp','wpwp','Inp','Dissip','Ret','yE','T','-v7.3')
if cspec==1
save([field_path,'statistics_sp_',fname,'.mat'],'Euu','Evv','Eww','Exx','Egg','Ezz','k','l','T')
end

%  save([field_path,'Umean_',fname,'.mat'],'Ust','Vst','Wst','Oxt','z','yE')
%  save([field_path,'statistics_spy_',fname,'.mat'],'Euu_y','Evv_y','Eww_y','Exx_y','Egg_y','Ezz_y','k','l','yE')
    


%%    


%load([field_path,'statistics_uu_',fname,'.mat'],'upvp','UP','upup','vpvp','wpwp','Inp','Dissip','Ret','y','T')
%load([field_path,'statistics_sp_',fname,'.mat'],'Euu','Evv','Eww','Exx','Egg','Ezz','k','l','T')
%rg=1:length(T); % range of fields used to calculate statistics

ys=1:(N+3)/2; % Half channel wall-normal grid
Retm=mean(Ret(1,rg)); % Average Re_{\tau}
ut=Retm/Re; % Average friction velocity

% Time-series
if cspec==1
E_comp % Farrell and Ioannou Energy plots.
M_comp % Hamilton, Kim and Waleffe RMS plots.
end

% One-point Statistics 
pt='og'; % Marker and Color 
if modf=='c'
    Stats_Couette
elseif modf=='m' || modf=='p'
    Stats_Poiseuille
end

