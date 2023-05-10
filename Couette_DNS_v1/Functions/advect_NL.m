function [adv_v,adv_g] = advect_NL(up,vp,wp,gp)


    dxdup=difX_F(up,1);
    dxdvp=difX_F(vp,1);
    dxdwp=difX_F(wp,1);
    %d2xdup=difX_F(up,2);
    %d2xdvp=difX_F(vp,2);
    %d2xdwp=difX_F(wp,2);
    
    
    dydup=difY_F(up,1);
    dydvp=difY_F(vp,1);
    dydwp=difY_F(wp,1);  
    %d2ydup=difY_F(up,2);
    %d2ydvp=difY_F(vp,2);
    %d2ydwp=difY_F(wp,2);
    
    
    dzdup=difZ_F(up,1);
    dzdvp=difZ_F(vp,1);
    dzdwp=difZ_F(wp,1);
    %d2zdup=difZ_F(up,2);
    %d2zdvp=difZ_F(vp,2);
    %d2zdwp=difZ_F(wp,2);
    
    %dxdgp=difX_F(gp,1);
    %dydgp=difY_F(gp,1);
    %dzdgp=difZ_F(gp,1);
    

    %D2gp=difX_F(gp,2)+difY_F(gp,2)+difZ_F(gp,2);

    %D2v=difX_F(vp,2)+difY_F(vp,2)+difZ_F(vp,2);

    %---------------------------------
    %v advection terms 
    ad_u=up.*dxdup+vp.*dydup+wp.*dzdup;
    ad_v=up.*dxdvp+vp.*dydvp+wp.*dzdvp;
    ad_w=up.*dxdwp+vp.*dydwp+wp.*dzdwp;
    %ad_g=up.*dxdgp+vp.*dydgp+wp.*dzdgp;
    
    %dzad_u=dzdup.*dxdup+dzdvp.*dydup+dzdwp.*dzdup;
    %dxad_w=dxdup.*dxdwp+dxdvp.*dydwp+dxdwp.*dzdwp;
    
    
    %ad1=difX_F(ad_v,2)+difZ_F(ad_v,2);%+difY_F(ad_v,2)
    ad1=-difX_F(ad_w,1)+difZ_F(ad_u,1);
    ad2=-difY_F(difX_F(ad_u,1)+difZ_F(ad_w,1),1)+(difX_F(ad_v,2)+difZ_F(ad_v,2));
   
    
    adv_v=ad2;

    
    adv_g=ad1;


end

