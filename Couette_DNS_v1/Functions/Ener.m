function [Etot] = Ener(F,G,H)
%Energy Density


global xE yE zE
%Extended Fields

%Etot=sum(F(:).^2+G(:).^2+H(:).^2)/NX/(N+1)/MZ/2

F(:,:,end+1)=F(:,:,1);
F(:,end+1,:)=F(:,1,:);

G(:,:,end+1)=G(:,:,1);
G(:,end+1,:)=G(:,1,:);

H(:,:,end+1)=H(:,:,1);
H(:,end+1,:)=H(:,1,:);

Etot=trapz(yE,trapz(xE,trapz(zE,F.*F+G.*G+H.*H,3),2),1)/2/xE(end)/2/zE(end);

% Euu=trapz(yE',trapz(xE,trapz(zE,F.*F,3),2),1)/2/xE(end)/yE(end)/2/zE(end);
% Evv=trapz(yE',trapz(xE,trapz(zE,G.*G,3),2),1)/2/xE(end)/yE(end)/2/zE(end);
% Eww=trapz(yE',trapz(xE,trapz(zE,H.*H,3),2),1)/2/xE(end)/yE(end)/2/zE(end);
% 
% Eout=[Euu;Evv;Eww];
end

