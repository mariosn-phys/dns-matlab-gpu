        %% Mean statistics Poiseuille
%        tst=101:4001;ys=1:37;
%        Retm=mean(Ret(tst));
%        ut=Retm/Re;
%     uu=squeeze(sum(sum(mean(Euu(2:end,1:end,:,rg),4),2),1)+sum(sum(mean(Euu(1,2:end,:,rg),4),2),1));
%     vv=squeeze(sum(sum(mean(Evv(2:end,1:end,:,rg),4),2),1)+sum(sum(mean(Evv(1,2:end,:,rg),4),2),1));
%     ww=squeeze(sum(sum(mean(Eww(2:end,1:end,:,rg),4),2),1)+sum(sum(mean(Eww(1,2:end,:,rg),4),2),1));
%     
%     uu=(uu+flip(uu))/2;
%     vv=(vv+flip(vv))/2;
%     ww=(ww+flip(ww))/2;
    
    upupf=mean(upup(:,rg),2);
    vpvpf=mean(vpvp(:,rg),2);
    wpwpf=mean(wpwp(:,rg),2);
    
    upupf=(upupf+flip(upupf))/2;
    vpvpf=(vpvpf+flip(vpvpf))/2;
    wpwpf=(wpwpf+flip(wpwpf))/2;
    
    
  %  figure;plot((y+1)*mean(Ret),[uu vv ww]*sqrt(2)/ut/ut) 
%     
%     for ij=1:7
%         
%         cc(:,ij)=coue171(:,ij);
%         
%     end
    load Re180.prof
    y2=Re180(:,2);
    U=Re180(:,3);
    u=Re180(:,4:6);
    uv=Re180(:,11);
    
    figure(321);clf
    semilogx(y2,U,'Linewidth',2)
    hold on
    plot((1+yE(ys))*Retm,mean(UP(ys,rg),2)/ut,pt)
    xlabel('$y^+$','Interpreter','latex')
    ylabel("$U /u_{\tau}$",'Interpreter','latex')
    
       figure(120);clf
       subplot(2,2,1)
    plot(y2,u(:,1),'Linewidth',2)
    hold on
    plot((1+yE(ys))*Retm,sqrt(mean(upupf(ys,:),2)/ut/ut),pt)
    %plot((1+yE(1:37))*Retm,2*sqrt(mean(uu(1:37,:),2)/ut/ut),'xk')
    xlabel('$y^+$','Interpreter','latex')
    ylabel("$u'/u_{\tau}$",'Interpreter','latex') 
    grid on
    
    
        
       %figure(120);clf
       subplot(2,2,2)
    plot(y2,u(:,2),'Linewidth',2)
    hold on
    plot((1+yE(ys))*Retm,sqrt(mean(vpvpf(ys,:),2)/ut/ut),pt)
    %plot((1+yE(1:37))*Retm,2*sqrt(mean(vv(1:37,:),2)/ut/ut),'xk')
    xlabel('$y^+$','Interpreter','latex')
    ylabel("$v'/u_{\tau}$",'Interpreter','latex') 
    grid on
    
        
    subplot(2,2,3)    
   %    figure(120);clf
    plot(y2,u(:,3),'Linewidth',2)
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
    plot((1+yE(ys))*Retm,mean(upvp(ys,rg),2)/ut/ut,pt)
    xlabel('$y^+$','Interpreter','latex')
    ylabel("$u'v'/u^2_{\tau}$",'Interpreter','latex')
    