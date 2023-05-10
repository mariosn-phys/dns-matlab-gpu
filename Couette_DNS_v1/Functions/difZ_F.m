function [dFdz] = difZ_F(F,n)
%%% Real fields
global D1z D2z N NX MZ

if n==1
    Dif=D1z;
elseif n==2
    Dif=D2z;
end

%dFdz=zeros(N+2,NX,MZ);

% for ix=1:NX
% dFdz(:,ix,:)=(Dif*permute(F(:,ix,:),[1 3 2])')';
% end

dFdz=permute(reshape(Dif*reshape(permute(F,[3 1 2]),[MZ,(N+2)*NX]),[MZ,N+2,NX]),[2 3 1]);

end

