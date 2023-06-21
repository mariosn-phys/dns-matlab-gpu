clear 
% New grid in y

addpath('../Functions/')

Nn=53; % Number of grid points in the wall-normal for the new flow field
% Fourier modes are adjusted when the flow field is loaded in the DNS code 
[~,yn]=cheb(Nn+1);yn=-yn;

resetT=0; % Set runtime of new file to 0

%Initial field name
field_path='../Data/Re600_n53/';Re=600;T=6148;
old_file=['state_Re',num2str(Re),'_',num2str(T(1),'%04.2f')];
save_new_path='../Data/';


load([field_path,old_file,'.mat'],'vhat','ghat','UP','WP','runtime')

Q=size(vhat); N=Q(1)-2;
[~,y]=cheb(N+1);
y=-y;

UP1n=interp1(y,UP,yn);
WP1n=interp1(y,WP,yn);

vhatn=zeros(Nn+2,Q(2),Q(3));
ghatn=vhatn;

for ik=1:Q(2)
    for ij=1:Q(3)
        
        vhatn(:,ik,ij)=interp1(y,vhat(:,ik,ij),yn);
        ghatn(:,ik,ij)=interp1(y,ghat(:,ik,ij),yn);

        %check interpolation          
        %figure(99);
        %subplot(1,2,1)
        %plot(imag(vhatn(:,ik,ij)),yn,'-*',imag(vhat(:,ik,ij)),y,'--')
        %subplot(1,2,2)
        %plot(imag(ghatn(:,ik,ij)),yn,'-*',imag(ghat(:,ik,ij)),y,'--')
        %drawnow
        %pause(.1)
    end
end

vhat=vhatn;
ghat=ghatn;
UP=UP1n;
WP=WP1n;

if resetT==1
    runtime=0;
    old_file=['state_Re',num2str(Re),'_',num2str(runtime,'%04.2f')];
end

save([save_new_path,old_file,'_new.mat'],'vhat','ghat','UP','WP','runtime')




