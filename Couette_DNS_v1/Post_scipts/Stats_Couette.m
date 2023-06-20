    %% Mean statistics

    
    upupf=mean(upup(:,rg),2);
    vpvpf=mean(vpvp(:,rg),2);
    wpwpf=mean(wpwp(:,rg),2);
    
    upupf=(upupf+flip(upupf))/2;
    vpvpf=(vpvpf+flip(vpvpf))/2;
    wpwpf=(wpwpf+flip(wpwpf))/2;
       
    
    c171=importdata('coue171.dat',' ',18);
    y2=c171.data(:,2);
    U=c171.data(:,3);
    u=c171.data(:,4);
    v=c171.data(:,5);
    w=c171.data(:,6);
    uv=c171.data(:,7);
    
    figure(321);clf
    semilogx(y2,U,'Linewidth',2)
    hold on
    plot((1+yE(ys))*Retm,mean(1+UP(ys,:),2)/ut,pt)
    xlabel('$y^+$','Interpreter','latex')
    ylabel("$U /u_{\tau}$",'Interpreter','latex')
    
       figure(120);clf
       subplot(2,2,1)
    plot(y2,u,'Linewidth',2)
    hold on
    plot((1+yE(ys))*Retm,sqrt(mean(upupf(ys,:),2)/ut/ut),pt)
    %plot((1+yE(1:37))*Retm,2*sqrt(mean(uu(1:37,:),2)/ut/ut),'xk')
    xlabel('$y^+$','Interpreter','latex')
    ylabel("$u'/u_{\tau}$",'Interpreter','latex') 
    grid on
    
    
        
       %figure(120);clf
       subplot(2,2,2)
    plot(y2,v,'Linewidth',2)
    hold on
    plot((1+yE(ys))*Retm,sqrt(mean(vpvpf(ys,:),2)/ut/ut),pt)
    %plot((1+yE(1:37))*Retm,2*sqrt(mean(vv(1:37,:),2)/ut/ut),'xk')
    xlabel('$y^+$','Interpreter','latex')
    ylabel("$v'/u_{\tau}$",'Interpreter','latex') 
    grid on
    
        
    subplot(2,2,3)    
   %    figure(120);clf
    plot(y2,w,'Linewidth',2)
    hold on
    plot((1+yE(ys))*Retm,sqrt(mean(wpwpf(ys,:),2)/ut/ut),pt)
    %plot((1+yE(1:37))*Retm,2*sqrt(mean(ww(1:37,:),2)/ut/ut),'xk')
    xlabel('$y^+$','Interpreter','latex')
    ylabel("$w'/u_{\tau}$",'Interpreter','latex') 
    grid on
    
    subplot(2,2,4)
    %figure(123);clf
    plot(y2,uv,'Linewidth',2)
    hold on
    plot((1+yE(ys))*Retm,mean(upvp(ys,:),2)/ut/ut,pt)
    xlabel('$y^+$','Interpreter','latex')
    ylabel("$u'v'/u^2_{\tau}$",'Interpreter','latex')
    