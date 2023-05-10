function [dFdy] = difY_F(F,n)

global DYF D2F N NX MZ

if n==1
    Dif=DYF;
elseif n==2
    Dif=D2F;
end

dFdy=zeros(N+2,NX,MZ);
 
% for iz=1:MZ
% dFdy(:,:,iz)=Dif*F(:,:,iz);
% end

dFdy=reshape(Dif*reshape(F,[N+2,NX*MZ]),[N+2,NX,MZ]);


end

