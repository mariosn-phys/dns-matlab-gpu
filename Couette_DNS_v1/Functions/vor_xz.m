function [gx,gz] = vor_xz(u0,v0,w0)
% x- and z- normal Vorticity fields 

% global a b N NX MZ D1x D2x D1z D2z DY D2 DYF D2F yE xE zE dy dx dz A B L 

dydw=difY_F(w0,1);
dzdv=difZ_F(v0,1);

dxdv=difX_F(v0,1);
dydu=difY_F(u0,1);

gx=dydw-dzdv;
gz=dxdv-dydu;


end

