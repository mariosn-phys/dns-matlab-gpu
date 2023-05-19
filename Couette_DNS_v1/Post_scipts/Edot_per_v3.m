function [E_MU_dot,E_s_dot,E_r_dot,E_p_dot,U_s_dot] = Edot_per_v3(u0,v0,w0)
%Energetics of elements in the mean - perturbation decomposition
%The terms 


global NX MZ yE zE L DYF D1z D2F D2z Re xE dy

    Um=repmat(mean(u0,2),[1,NX,1]);Us=squeeze(mean(Um-repmat(mean(mean(Um,3),2),[1,NX,MZ]),2));UP=squeeze(mean(mean(u0,3),2));
    upert=u0-Um;
    Vm=repmat(mean(v0,2),[1,NX,1]);Vs=squeeze(mean(Vm-repmat(mean(mean(Vm,3),2),[1,NX,MZ]),2));
    vpert=v0-Vm;
    Wm=repmat(mean(w0,2),[1,NX,1]);Ws=squeeze(mean(Wm-repmat(mean(mean(Wm,3),2),[1,NX,MZ]),2));WP=squeeze(mean(mean(w0,3),2));
    wpert=w0-Wm;
    
    
    dxdup=difX_F(upert,1);
    dxdvp=difX_F(vpert,1);
    dxdwp=difX_F(wpert,1);
       
    dydup=difY_F(upert,1);
    dydvp=difY_F(vpert,1);
    dydwp=difY_F(wpert,1);  
       
    dzdup=difZ_F(upert,1);
    dzdvp=difZ_F(vpert,1);
    dzdwp=difZ_F(wpert,1);
    
    ad_u=upert.*dxdup+vpert.*dydup+wpert.*dzdup;
    ad_v=upert.*dxdvp+vpert.*dydvp+wpert.*dzdvp;
    ad_w=upert.*dxdwp+vpert.*dydwp+wpert.*dzdwp;
    
    dxdU=difX_F(Um,1);
    dxdV=difX_F(Vm,1);
    dxdW=difX_F(Wm,1);
       
    dydU=difY_F(Um,1);
    dydV=difY_F(Vm,1);
    dydW=difY_F(Wm,1);  
       
    dzdU=difZ_F(Um,1);
    dzdV=difZ_F(Vm,1);
    dzdW=difZ_F(Wm,1);
    
    ad_U=Um.*dxdU+Vm.*dydU+Wm.*dzdU;
    ad_V=Um.*dxdV+Vm.*dydV+Wm.*dzdV;
    ad_W=Um.*dxdW+Vm.*dydW+Wm.*dzdW;
    
    %stresses
    UmVm=squeeze(mean(Um.*Vm,2));
    UmVm=UmVm-repmat(mean(UmVm,2),[1,MZ]);
    UsVs=Us.*Vs;UsVs_p=squeeze(mean(UsVs,2));
    UsVs=UsVs-repmat(mean(UsVs,2),[1,MZ]);
    
    UmWm=squeeze(mean(Um.*Wm,2));
    UmWm=UmWm-repmat(mean(UmWm,2),[1,MZ]);
    UsWs=Us.*Ws;
    UsWs=UsWs-repmat(mean(UsWs,2),[1,MZ]);
    
    u_pv_p=squeeze(mean(upert.*vpert,2));u_pv_p_p=squeeze(mean(u_pv_p,2));
    u_pv_p=u_pv_p-repmat(mean(u_pv_p,2),[1,MZ]);    
    u_pw_p=squeeze(mean(upert.*wpert,2));
    u_pw_p=u_pw_p-repmat(mean(u_pw_p,2),[1,MZ]);
    v_pv_p=squeeze(mean(vpert.*vpert,2));
    
    w_pw_p=squeeze(mean(wpert.*wpert,2));
    
    v_pw_p=squeeze(mean(vpert.*wpert,2));v_pw_p_p=squeeze(mean(v_pw_p,2));
    
    
    %diff
    dUmVmdY=DYF*UmVm;
    dUsVsdY=DYF*UsVs;
    dUmWmdZ=(D1z*(UmWm)')';
    dUsWsdZ=(D1z*(UsWs)')';
    dupvpdY=DYF*u_pv_p;
    dupwpdZ=(D1z*(u_pw_p)')';
    D2Us=D2F*Us+Us*D2z;
    dUsdY=DYF*Us;
    dUsdZ=(D1z*(Us)')';
    
    dVmdY=DYF*Vs;
    dVmdZ=(D1z*(Vs)')';
    dWmdY=DYF*Ws;
    D2Vs=D2F*Vs+Vs*D2z;
    D2Ws=D2F*Ws+Ws*D2z;
    
    
    dvpvpdY=DYF*v_pv_p;
    dvpwpdZ=(D1z*(v_pw_p)')';
    dvpwpdY=DYF*v_pw_p;
    dwpwpdZ=(D1z*(w_pw_p)')';
    D2u=difX_F(upert,2)+difY_F(upert,2)+difZ_F(upert,2);
    D2v=difX_F(vpert,2)+difY_F(vpert,2)+difZ_F(vpert,2);
    D2w=difX_F(wpert,2)+difY_F(wpert,2)+difZ_F(wpert,2);
    
    dyDU=DYF*UP;%dyDU(1)=(UP(2)-UP(1))/dy;dyDU(end)=(UP(end)-UP(end-1))/dy;
    dyDW=DYF*WP;%dyDW(1)=(WP(2)-WP(1))/dy;dyDW(end)=(WP(end)-WP(end-1))/dy;
       
    
    %last point in z 
    Us=[Us Us(:,1)];    
    dUmVmdY=[dUmVmdY dUmVmdY(:,1)];
    dUsVsdY=[dUsVsdY dUsVsdY(:,1)];
    dUmWmdZ=[dUmWmdZ dUmWmdZ(:,1)];
    dUsWsdZ=[dUsWsdZ dUsWsdZ(:,1)];
    
    dupvpdY=[dupvpdY dupvpdY(:,1)];
    dupwpdZ=[dupwpdZ dupwpdZ(:,1)];
    D2Us=[D2Us D2Us(:,1)];
    dUsdY=[dUsdY dUsdY(:,1)];
    dUsdZ=[dUsdZ dUsdZ(:,1)];
    UmVm=[UmVm UmVm(:,1)];
    UmWm=[UmWm UmWm(:,1)];
    
    Vs=[Vs Vs(:,1)];
    Ws=[Ws Ws(:,1)];
    
    
    v_pv_p=[v_pv_p v_pv_p(:,1)];  
    w_pw_p=[w_pw_p w_pw_p(:,1)]; 
    v_pw_p=[v_pw_p v_pw_p(:,1)];
    dVmdY=[dVmdY dVmdY(:,1)];
    dVmdZ=[dVmdZ dVmdZ(:,1)];
    dWmdY=[dWmdY dWmdY(:,1)];
    D2Vs=[D2Vs D2Vs(:,1)];
    D2Ws=[D2Ws D2Ws(:,1)];
    
    dvpvpdY=[dvpvpdY dvpvpdY(:,1)];
    dvpwpdZ=[dvpwpdZ dvpwpdZ(:,1)];
    dvpwpdY=[dvpwpdY dvpwpdY(:,1)];
    dwpwpdZ=[dwpwpdZ dwpwpdZ(:,1)];
    
    Um(:,:,end+1)=Um(:,:,1);
    Um(:,end+1,:)=Um(:,1,:);
    upert(:,:,end+1)=upert(:,:,1);
    upert(:,end+1,:)=upert(:,1,:);    
    D2u(:,:,end+1)=D2u(:,:,1);
    D2u(:,end+1,:)=D2u(:,1,:);

    vpert(:,:,end+1)=vpert(:,:,1);
    vpert(:,end+1,:)=vpert(:,1,:); 
    D2v(:,:,end+1)=D2v(:,:,1);
    D2v(:,end+1,:)=D2v(:,1,:);

    wpert(:,:,end+1)=wpert(:,:,1);
    wpert(:,end+1,:)=wpert(:,1,:); 
    D2w(:,:,end+1)=D2w(:,:,1);
    D2w(:,end+1,:)=D2w(:,1,:);
   
    ad_u(:,:,end+1) = ad_u(:,:,1);
    ad_u(:,end+1,:) = ad_u(:,1,:); 
    ad_v(:,:,end+1) = ad_v(:,:,1);
    ad_v(:,end+1,:) = ad_v(:,1,:);
    ad_w(:,:,end+1) = ad_w(:,:,1);
    ad_w(:,end+1,:) = ad_w(:,1,:);
    
    ad_U(:,:,end+1) = ad_U(:,:,1);
    ad_U(:,end+1,:) = ad_U(:,1,:); 
    ad_V(:,:,end+1) = ad_V(:,:,1);
    ad_V(:,end+1,:) = ad_V(:,1,:);
    ad_W(:,:,end+1) = ad_W(:,:,1);
    ad_W(:,end+1,:) = ad_W(:,1,:);
    
%     %Terms
%     
    E_s_dot(1)=trapz(yE,trapz(zE,Us.*D2Us,2))/Re/L/zE(end);
    E_s_dot(2)=trapz(yE,trapz(zE,-Us.*Vs.*(repmat(DYF*UP,[1 MZ+1])),2))/L/zE(end);
    E_s_dot(3)=trapz(yE,trapz(zE,-Us.*dUsVsdY,2))/L/zE(end);
    E_s_dot(4)=trapz(yE,trapz(zE,-Us.*dUsWsdZ,2))/L/zE(end);
    E_s_dot(5)=trapz(yE,trapz(zE,-Us.*dupvpdY,2))/L/zE(end);
    E_s_dot(6)=trapz(yE,trapz(zE,-Us.*dupwpdZ,2))/L/zE(end);

%     
    U_s_dot(1)=trapz(yE,trapz(zE,sign(Us).*D2Us,2))/Re/L/zE(end);
    U_s_dot(2)=trapz(yE,trapz(zE,-sign(Us).*Vs.*(repmat(DYF*UP,[1 MZ+1])),2))/L/zE(end);
    U_s_dot(3)=trapz(yE,trapz(zE,-sign(Us).*dUsVsdY,2))/L/zE(end);
    U_s_dot(4)=trapz(yE,trapz(zE,-sign(Us).*dUsWsdZ,2))/L/zE(end);
    U_s_dot(5)=trapz(yE,trapz(zE,-sign(Us).*dupvpdY,2))/L/zE(end);
    U_s_dot(6)=trapz(yE,trapz(zE,-sign(Us).*dupwpdZ,2))/L/zE(end);
%     
%     E_s_dot_a(1)=trapz(yE,trapz(zE,dUsdZ.*UmWm,2))/L/zE(end);
%     E_s_dot_a(2)=trapz(yE,trapz(zE,-Us.*Vs.*(repmat(DYF*UP,[1 MZ+1])),2))/L/zE(end);
%     
%     Max_Usz=max(max(dUsdZ));
%     Max_Usy=max(max(dUsdY));
%     
%     
    E_r_dot(1)=trapz(yE,trapz(zE,Vs.*D2Vs+Ws.*D2Ws,2))/Re/L/zE(end);
    E_r_dot(2)=trapz(yE,trapz(zE,(v_pv_p-w_pw_p).*dVmdY+v_pw_p.*(dVmdZ+dWmdY),2))/L/zE(end);
%     %E_r_dot(3)=trapz(yE,trapz(zE,v_pw_p.*(dVmdZ+dWmdY),2))/L/zE(end);
%     
    E_p_dot(1)=trapz(yE,trapz(xE,trapz(zE,upert.*D2u+vpert.*D2v+wpert.*D2w,3),2))/Re/L/zE(end)/xE(end);
    E_p_dot(2)=trapz(yE,trapz(zE,Vs.*dvpvpdY+Ws.*dvpwpdY,2))/L/zE(end);
    E_p_dot(3)=trapz(yE,trapz(zE,Vs.*dvpwpdZ+Ws.*dwpwpdZ,2))/L/zE(end);
    E_p_dot(4)=trapz(yE,UP.*(DYF*u_pv_p_p))/L;
    E_p_dot(5)=trapz(yE,WP.*(DYF*v_pw_p_p))/L;
    %E_p_dot(3)=trapz(yE,trapz(zE,Ws.*dvpwpdY+Ws.*dwpwpdZ,2))/L/zE(end);
    E_p_dot(6)=trapz(yE,trapz(zE,Us.*dupvpdY,2))/L/zE(end);
    E_p_dot(7)=trapz(yE,trapz(zE,Us.*dupwpdZ,2))/L/zE(end);
    E_p_dot(8)=-trapz(yE,trapz(xE,trapz(zE,upert.*ad_u+vpert.*ad_v+wpert.*ad_w,3),2))/L/zE(end)/xE(end);
%        
%     %E_MU_dot(1)=(dyDU(1)+dyDU(end))/Re/L-trapz(yE,dyDU.^2)/Re/L;
     E_MU_dot(1)=-trapz(yE,dyDU.^2)/Re/L;
     E_MU_dot(2)=(dyDU(1)+dyDU(end))/Re/L;
     E_MU_dot(3)=-trapz(yE,trapz(xE,trapz(zE,Um.*ad_U,3),2))/L/zE(end)/xE(end);
     E_MU_dot(4)=-trapz(yE,trapz(xE,trapz(zE,Um.*ad_u,3),2))/L/zE(end)/xE(end);

%     
%     %E_MW_dot(1)=trapz(yE,dyDW.^2)/Re/L;
%     %E_MW_dot(2)=trapz(yE,dyDW(1)+dyDW(end))/Re/L;
    
    
    
end

