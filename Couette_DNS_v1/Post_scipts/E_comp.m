    %% 

    
    % Mean-Streak-Roll-perturbation(kx) Variables
    EM=[];
    EUs=[];
    Er=[];
    Ep=[];
    Ekx=[];
    
    %rg=1:size(Euu,3);Nitem=length(rg);Thalf=0:floor(Nitem/2);
    tt=T(rg);
    EM=[EM;squeeze(Euu(1,1,rg))];
    EUs=[EUs;squeeze(sum(Euu(1,2:end,rg),2))];
    Er=[Er;squeeze(sum(Evv(1,2:end,rg)+Eww(1,2:end,rg),2))];
    Ep=[Ep;squeeze(sum(sum(Euu(2:end,:,rg)+Evv(2:end,:,rg)+Eww(2:end,:,rg),1),2))];
    Ekx=[Ekx;squeeze(sum(Euu(2:end,:,rg)+Evv(2:end,:,rg)+Eww(2:end,:,rg),2))];
    
    figure(111);clf
    plot(tt,Ekx,'LineWidth',2)
    xlabel('$t U_w/h$','Interpreter','Latex')
    ylabel('$E_{k_x}$','Interpreter','Latex')
    set(gca,'LineWidth',2,'FontSize',16)
    %xlim([3500 4000])
    
    figure(112);clf
    plot(k(2:size(Ekx)+1),mean(Ekx,2),'-+','Linewidth',2)
    xlabel('$k_x$','Interpreter','Latex')
    ylabel('$<E_{k_x}>$','Interpreter','Latex')
    set(gca,'LineWidth',2,'FontSize',16)
   
    figure(16);
    col=[0 0 1;1 0 0;0.16 0.5 0.15;0 0 0];
    box on
    set(gca, 'ColorOrder', col, 'NextPlot', 'replacechildren','LineWidth',2);
    %h2=semilogy(T,[EM EUs Er Ep Ex Ey Ez]);
    h2=semilogy(tt,[EM EUs Er Ep]);

    set(h2,'LineWidth',2)
    xlabel('$t U_w/h$','Interpreter','Latex')
    ylabel('$E_i$','Interpreter','Latex')
    
    h1=legend("$E_M$","$E_U$","$E_r$","$E_p$",'Orientation','horizontal');
    set(h1,'Interpreter','latex')
    
    col1=[0 0.4470 0.7410;0.8500    0.3250    0.0980;0.9290    0.6940    0.1250;0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880;0.3010    0.5950    0.9330;0.1 0.9 0.1;0 0 0 ];

%    save([field_path,'Energies_',fname,'.mat'],'EM','EUs','Er','Ep','Ekx','T')
    
    