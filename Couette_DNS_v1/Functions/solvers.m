function [ICvkron,ICvDvkron,ICgkron,ICggkron,ICvCv,ICvIDELv,ICg,ICgCg,S_mf,S_mp,Sol_m,kkm,llm]=solvers(g,h)

%%% Mean flow and v,g solver matrices with timestep of fraction h

global N NX MZ D2 D2F k l D4

    g=g*h;

    %Precalc Matrices I and Poisson Solvers
    I=eye(N);
    I2=eye(N+2);
    
    %SolU
    %SolW
    
    S_mm=I-g*D2;
    S_mf=I2-g*D2F;
    S_mp=I2+g*D2F;
    Sol_m=I/S_mm;
    
    
    ICg=zeros(N,N,NX,MZ);
    Cg=zeros(N,N,NX,MZ);
    
    ICv=zeros(N,N,NX,MZ);
    Cv=zeros(N,N,NX,MZ);
    
    
  for ii=1:NX/3
    kp=k(ii);
    for kk=1:length(l)
    lp=l(kk);
     
    alpha=kp^2+lp^2;    
    
    DEL=D2-alpha*I;
    IDEL=DEL\I;
    	
	IDELv(:,:,ii,kk)=IDEL;
    
    Scn_m=I-g*DEL;
    Scn_p=I+g*DEL;
        
    ICg(:,:,ii,kk)=Scn_m\I;
    %Cg(:,:,ii,kk)=Scn_p;
        
    ICgCg(:,:,ii,kk)=ICg(:,:,ii,kk)*Scn_p;
    
    DEL4=D4-2*alpha*D2+2*kp^2*lp^2*I+kp^4*I+lp^4*I;
    
    Svim=I-g*IDEL*DEL4;
    Svip=I+g*IDEL*DEL4;

    ICv=Svim\I;
    %ICv(:,:,ii,kk)=Svim\I;
    %Cv(:,:,ii,kk)=Svip;
    
    ICvCv(:,:,ii,kk)=ICv*Svip;
    ICvIDELv(:,:,ii,kk)=ICv*IDEL;
    end
  end
  
          ICvkron=spalloc(N*(NX/3 )*(2*MZ/3-1),N*(NX/3 )*(2*MZ/3-1),N*N*(NX/3 )*(2*MZ/3-1));
          ICvDvkron=ICvkron;
          ICgkron=ICvkron;
          ICggkron=ICvkron;
          kkm=spalloc(N*(NX/3 )*(2*MZ/3-1),N*(NX/3 )*(2*MZ/3-1),N*(NX/3 )*(2*MZ/3-1));
          llm=kkm;

% MZ,NX
for ii=1:NX/3 
    tic;
    
    ICvkronb=ICvCv(:,:,ii,1);
    ICvDvkronb=ICvIDELv(:,:,ii,1);
    ICgkronb=ICg(:,:,ii,1);
    ICggkronb=ICgCg(:,:,ii,1);
    
    kpv=k(ii);
    lpv=l(1);
    alp=kpv^2+lpv^2;
    
    lpb=zeros(N);
    if ii==1
        kpb=zeros(N);
    else
        kpb=1i*kpv*eye(N)/alp;
    end

    
    for jk=[2:MZ/3  MZ-MZ/3+2:MZ]
        
        ICvkronb=blkdiag(ICvkronb,ICvCv(:,:,ii,jk));
        ICvDvkronb=blkdiag(ICvDvkronb,ICvIDELv(:,:,ii,jk));
        ICgkronb=blkdiag(ICgkronb,ICg(:,:,ii,jk));
        ICggkronb=blkdiag(ICggkronb,ICgCg(:,:,ii,jk));
        
        lpv=l(jk);
        alp=kpv^2+lpv^2;
        kpb=blkdiag(kpb,1i*kpv*eye(N)/alp);
        lpb=blkdiag(lpb,1i*lpv*eye(N)/alp);
        
    end
    
    ei=zeros(NX/3,1);ei(ii)=1;
    
    IM=spdiags(ei,0,NX/3,NX/3);
    ICvkron=ICvkron+kron(IM,ICvkronb);
    ICvDvkron=ICvDvkron+kron(IM,ICvDvkronb);
    ICgkron=ICgkron+kron(IM,ICgkronb);
    ICggkron=ICggkron+kron(IM,ICggkronb);
    
    kkm=kkm+kron(IM,kpb);
    llm=llm+kron(IM,lpb);
    toc
end
      
  
%  NX,MZ           
%  jit=0;
%   for jk=[1:MZ/3+1 MZ-MZ/3+1:MZ]
%       jit=jit+1;
%       
%       ICvkronb=ICvCv(:,:,1,jk);
%       ICvDvkronb=ICvIDELv(:,:,1,jk);
%       ICgkronb=ICg(:,:,1,jk);
%       ICggkronb=ICgCg(:,:,1,jk);
%       
%       for ii=2:NX/3+1
%           ICvkronb=blkdiag(ICvkronb,ICvCv(:,:,ii,jk));
%           ICvDvkronb=blkdiag(ICvDvkronb,ICvIDELv(:,:,ii,jk));
%           ICgkronb=blkdiag(ICgkronb,ICg(:,:,ii,jk));
%           ICggkronb=blkdiag(ICggkronb,ICgCg(:,:,ii,jk));  
%       end
%       
%       ei=zeros(2*MZ/3+1,1);ei(jit)=1;
%       
%       IM=spdiags(ei,0,2*MZ/3+1,2*MZ/3+1);
%       ICvkron=ICvkron+kron(IM,ICvkronb);
%       ICvDvkron=ICvDvkron+kron(IM,ICvDvkronb);
%       ICgkron=ICgkron+kron(IM,ICgkronb);
%       ICggkron=ICggkron+kron(IM,ICggkronb);
%       
%   end
%   
  
  
  
  
  
  
  
  
  
