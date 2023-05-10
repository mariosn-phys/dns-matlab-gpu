clear all;

global M N Re NX D1x D2x D1z D2z DYF D2F D2 D4 y gamma
global MZ k l dt xE yE zE A B L dx dz
global b1 b2
global kLyap NLyap
global a b kkm llm kkmm llmm DYkron DYkronm
global x yE z
global xx_yx yy_yx yy_yz zz_yz xx_xz zz_xz

%gpuDevice(1)


addpath('../Functions/')
%addpath('Functions_metrics/')
%addpath('Functions_modal/')
%field_path='../DNS_cheb_RK3_v3/Re600_cheb_modal_uvw_n53_rnl_control3_dconv/'

field_path='../Data/Re600_n53/'  %% save and restart file path


Re=600;%Re=dum(2);
% Poiseuille or Couette 'p' or 'c'
mod='c';
it0=2;
dt=0.02;
Ti=0;Tf=6155;Tstp=1;
T=Ti:Tstp:Tf;
NT=length(T);
g=dt/(2*Re);

a=2/(1.75);%a=2*pi/1.774;
b=2/(1.2);%b=2*pi/0.884;

%%%% NX,MZ 12,24,48, +12 or +6

A=2*pi/a;NX=54;dx=A/NX;x=-dx+dx*(1:NX)';xE=[x;A]; 
B=2*pi/b;MZ=54;dz=B/MZ;z=-dz+dz*(1:MZ)';zE=[z;B];
L=2;

N=53;


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
  
    if mod=='p'
    U1=(1-yE.^2);gm=-2/Re; %Poiseuille
    elseif mod=='c'
    U1=yE;gm=0; %Couette
    elseif mod=='z'
    U1=yE*0;gm=-2/Re*0;    
    end
    b1=U1(1);b2=U1(end);

    %Initial Conditions vortex+optimal n_x=1,n_z=1
  
    Uback=repmat(U1,[1,NX,MZ]);
    gamma=gm*ones(N+2,1);

  
   iter=0;
   
   yp=21; % y cross section to plot

    for it=1:NT
    iter=iter+1
    
    start_file=[field_path,'state_Re',num2str(Re),'_',num2str(T(it),'%04.2f'),'.mat']

   [vi,gi,UP1,WP1]=read_from_disk_compact(start_file);
   [u0,v0,w0,g0] = make_uw(gather(vi),gather(gi),gather(UP1),gather(WP1),b1,b2);
 
    measure_plot_ener
    upvp(:,iter)=mean(mean(u0.*v0,3),2);
    UP(:,iter)=gather(UP1);
   

    plot_mode_2(17,u0,v0,w0,g0,22,12,16) % Plot velocity state
    drawnow
    
%     vmean1=squeeze(mean(v0(:,:,:,1),2)); % vmean=0*vmean;
%     gmean1=squeeze(mean(g0(:,:,:,1),2)); % gmean=0*gmean;
%     umean1=squeeze(mean(u0(:,:,:,1),2));
%     wmean1=squeeze(mean(w0(:,:,:,1),2));
%     
%     Ust(:,:,iter)=umean1;
%     Vst(:,:,iter)=vmean1;
%     Wst(:,:,iter)=wmean1;
%     Oxt(:,:,iter)=squeeze(mean(difY_F(w0,1)-difZ_F(v0,1),2));
       
    end
    
   fname=['R600'];

%    EM=[];
%    EUs=[];
%    Er=[];
%    Ep=[];Ekx=[];
%    
%    rg=1:iter;
%    


%  save([field_path,'Umean_',fname,'.mat'],'Ust','Vst','Wst','Oxt','z','yE')
  save([field_path,'statistics_',fname,'.mat'],'upvp','UP','Euu','Evv','Eww','Egg','Exx','Ezz','Efm','Dissip','Inp','O_bot','O_top','T')

    
    %% 

    
    % Mean-Streak-Roll-perturbation(kx) Variables
    EM=[];
    EUs=[];
    Er=[];
    Ep=[];
    Ekx=[];
    
    rg=1:size(Euu,3);Nitem=length(rg);Thalf=0:floor(Nitem/2);
   
    EM=[EM;squeeze(Euu(1,1,rg))];
    EUs=[EUs;squeeze(sum(Euu(1,2:end,rg),2))];
    Er=[Er;squeeze(sum(Evv(1,2:end,rg)+Eww(1,2:end,rg),2))];
    Ep=[Ep;squeeze(sum(sum(Euu(2:end,:,rg)+Evv(2:end,:,rg)+Eww(2:end,:,rg),1),2))];
    Ekx=[Ekx;squeeze(sum(Euu(2:end,:,rg)+Evv(2:end,:,rg)+Eww(2:end,:,rg),2))];
    
    figure(111);clf
    plot(T,Ekx,'LineWidth',2)
    xlabel('$t U_w/h$','Interpreter','Latex')
    ylabel('$E_{k_x}$','Interpreter','Latex')
    set(gca,'LineWidth',2,'FontSize',16)
    %xlim([3500 4000])
    
    figure(112);clf
    plot(k(2:size(Ekx)+1),mean(Ekx,2),'-+','Linewidth',2)
    xlabel('$k_x$','Interpreter','Latex')
    ylabel('$<E_{k_x}>$','Interpreter','Latex')
    set(gca,'LineWidth',2,'FontSize',16)

    
    Ret=sqrt(Inp*Re);
    ut=mean(Ret)/Re;
    
   
    figure(16);
    col=[0 0 1;1 0 0;0.16 0.5 0.15;0 0 0];
    box on
    set(gca, 'ColorOrder', col, 'NextPlot', 'replacechildren','LineWidth',2);
    %h2=semilogy(T,[EM EUs Er Ep Ex Ey Ez]);
    h2=semilogy(T,[EM EUs Er Ep]);

    set(h2,'LineWidth',2)
    xlabel('$t U_w/h$','Interpreter','Latex')
    ylabel('$E_i$','Interpreter','Latex')
    
    h1=legend("$E_M$","$E_U$","$E_r$","$E_p$",'Orientation','horizontal');
    set(h1,'Interpreter','latex')
    
    col1=[0 0.4470 0.7410;0.8500    0.3250    0.0980;0.9290    0.6940    0.1250;0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880;0.3010    0.5950    0.9330;0.1 0.9 0.1;0 0 0 ];

%    save([field_path,'Energies_',fname,'.mat'],'EM','EUs','Er','Ep','Ekx','T')
    
    
    %% Hamilton Kim Waleffe
    % Square root of energy at selected wavenumber pair
    
    M0b=[];
    M02b=[];
    M03b=[];
    Ma0=[];
    M2a0=[];
    Mab=[];
    Ma2b=[];
    Ma3b=[];
    M2ab=[];
    M2a2b=[];
        
    rg=1:size(Euu,3);Nitem=length(rg);Thalf=0:floor(Nitem/2);
    M0b=[M0b;squeeze(Euu(1,2,rg)+Evv(1,2,rg)+Eww(1,2,rg))];
    M02b=[M02b;squeeze(Euu(1,3,rg)+Evv(1,3,rg)+Eww(1,3,rg))];
    M03b=[M03b;squeeze(Euu(1,4,rg)+Evv(1,4,rg)+Eww(1,4,rg))];
    Ma0=[Ma0;squeeze(Euu(2,1,rg)+Evv(2,1,rg)+Eww(2,1,rg))];
    M2a0=[M2a0;squeeze(Euu(3,1,rg)+Evv(3,1,rg)+Eww(3,1,rg))];    
    Mab=[Mab;squeeze(Euu(2,2,rg)+Evv(2,2,rg)+Eww(2,2,rg))];
    Ma2b=[Ma2b;squeeze(Euu(2,3,rg)+Evv(2,3,rg)+Eww(2,3,rg))];
    Ma3b=[Ma3b;squeeze(Euu(2,4,rg)+Evv(2,4,rg)+Eww(2,4,rg))];
    M2ab=[M2ab;squeeze(Euu(3,2,rg)+Evv(3,2,rg)+Eww(3,2,rg))];
    M2a2b=[M2a2b;squeeze(Euu(3,3,rg)+Evv(3,3,rg)+Eww(3,3,rg))];
  

    figure(12);semilogy(T,sqrt([M0b M02b M03b Ma0 M2a0 Mab Ma2b Ma3b M2ab M2a2b]),'LineWidth',2)
    hl=legend("$M_{(0,b)}$","$M_{(0,2b)}$","$M_{(0,3b)}$","$M_{(a,0)}$","$M_{(2a,0)}$","$M_{(a,b)}$","$M_{(a,2b)}$","$M_{(a,3b)}$","$M_{(2a,b)}$","$M_{(2a,2b)}$",'Orientation','horizontal');
    set(hl,'Interpreter','latex')
    xlabel('$t$','Interpreter','Latex')
    ylabel('$M_{(a,b)}$','Interpreter','Latex')
    set(gca,'LineWidth',2,'FontSize',16)
    
