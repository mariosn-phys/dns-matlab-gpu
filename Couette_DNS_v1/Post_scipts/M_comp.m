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
        
    tt=T(rg);
    %rg=1:size(Euu,3);Nitem=length(rg);Thalf=0:floor(Nitem/2);
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
  

    figure(12);semilogy(tt,sqrt([M0b M02b M03b Ma0 M2a0 Mab Ma2b Ma3b M2ab M2a2b]),'LineWidth',2)
    hl=legend("$M_{(0,b)}$","$M_{(0,2b)}$","$M_{(0,3b)}$","$M_{(a,0)}$","$M_{(2a,0)}$","$M_{(a,b)}$","$M_{(a,2b)}$","$M_{(a,3b)}$","$M_{(2a,b)}$","$M_{(2a,2b)}$",'Orientation','horizontal');
    set(hl,'Interpreter','latex')
    xlabel('$t$','Interpreter','Latex')
    ylabel('$M_{(a,b)}$','Interpreter','Latex')
    set(gca,'LineWidth',2,'FontSize',16)