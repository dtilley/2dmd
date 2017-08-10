% initalize input variables 
    n=20; % square lattice edge length
    N=n*n; % total lattice sites
    t=1000000; % number of MC steps
    in=7/8; % fraction of inital atoms with downspin
    dspin=round(N*in);
    H=0; % external magnetic field 
    T=0.5:0.1:2.8; % temperature vector
    mu=1; % magnetic moment
    J=1; % coupling coefficient
    prefix='isingMC-20x20-J1H0'; % output file prefix
    wrt=100; % write out every wrt steps
    equil=0.35; % throw out the first 35% of the data
% initalize output variables
    s=zeros(1,length(T));
    Cv=zeros(1,length(T));
for k=1:length(T)
% in:the starting condition in fraction of downspins
    energy=zeros(1,round(t/wrt));
    mag=zeros(1,round(t/wrt));
    E=0;
    % M=0;
    l=ones(n,n);
    xdspin=round(rand(1,dspin)*n);
    ydspin=round(rand(1,dspin)*n);
    isunique=0;
    while(~isunique)
        for i=1:length(xdspin)
           if(xdspin(i)==0) % check if the x-coord is zero
            xdspin(i)=round(rand(1)*n);
           end
           if(ydspin(i)==0) % check if the y-coord is zero
            ydspin(i)=round(rand(1)*n);
           end
           % check if xy-coord is unique
           count=0;
           for j=i+1:length(xdspin)
            if(xdspin(i)==xdspin(j) && ydspin(i)==ydspin(j))
                xdspin(j)=round(rand(1)*n);
                ydspin(j)=round(rand(1)*n);
                count=count+1;
            end
           end
        end
        if (count==0 && isempty(find(xdspin==0)) && isempty(find(ydspin==0)))
            isunique=1;
        end
    end
    
    for i=1:length(xdspin)
       l(xdspin(i),ydspin(i))=-1; 
    end
% Calculate E & M initial
for i=1:n
   for j=1:n
      if (i==1 && j~=n) % non-corner top edge
          E=E+-H*mu*l(i,j)-J*(l(i,j)*l(n,j)+l(i,j)*l(i,j+1));
      elseif (i==1 && j==n) % upper right corner
          E=E+-H*mu*l(i,j)-J*(l(i,j)*l(n,j)+l(i,j)*l(i,1));
      elseif (i~=1 && j==n) % non-corner right edge
          E=E+-H*mu*l(i,j)-J*(l(i,j)*l(i,1)+l(i,j)*l(i-1,j));
      else % non-boundary neighbor above and to the right
          E=E+-H*mu*l(i,j)-J*(l(i,j)*l(i-1,j)+l(i,j)*l(i,j+1));
      end 
   end
end
M=abs(sum(sum(l)));
% Monte Carlo steps
ndx=1;
for i=1:t
    if(mod(i,wrt)==0)
        energy(ndx)=E;
        mag(ndx)=abs(M);
        ndx=ndx+1;
    end
    isnonzero=0;
    xMC=round(rand(1)*n);
    yMC=round(rand(1)*n);
    while (~isnonzero)
        if (xMC~=0 && yMC~=0)
           isnonzero=1;
        elseif (xMC==0)
            xMC=round(rand(1)*n);
        else
            yMC=round(rand(1)*n);
        end
    end
   % calculate deltaE
   sMC=l(xMC,yMC)*-1; % flip the spin
   deltaE=-2*H*mu*sMC; % calculate external field energy
   dEspin=0;
   % calculate spin couple energy row coord
   if (xMC==n) 
    dEspin=dEspin+(l(1,yMC)+l(xMC-1,yMC));
   elseif (xMC==1)
    dEspin=dEspin+(l(xMC+1,yMC)+l(n,yMC));
   else
    dEspin=dEspin+(l(xMC+1,yMC)+l(xMC-1,yMC));
   end
   % calculate spin couple energy col coord
   if (yMC==n)
    dEspin=dEspin+(l(xMC,1)+l(xMC,yMC-1));
   elseif (yMC==1)
    dEspin=dEspin+(l(xMC,n)+l(xMC,yMC+1));
   else
    dEspin=dEspin+(l(xMC,yMC+1)+l(xMC,yMC-1));
   end
   deltaE=deltaE-J*sMC*dEspin;
   if (deltaE<0) % if the change if favorable accept
       l(xMC,yMC)=sMC;
       %M=M+2*mu*sMC;
       M=abs(sum(sum(l)));
       E=E+deltaE;
   elseif (exp(-deltaE/T(k))>=rand(1)) % if the probablity of finding that energy
       l(xMC,yMC)=sMC;
       %M=M+2*mu*sMC;
       M=abs(sum(sum(l)));
       E=E+deltaE;
   end
end

% calculate Cv and <m>
mag=mag(round(equil*length(mag)):end);
energy=energy(round(equil*length(energy)):end);
meanEsqr=mean(energy.*energy);
meanE=mean(energy);

s(k)=mean(mag)/N;
Cv(k)=(meanEsqr-meanE*meanE)/(N*T(k)*T(k));

% write final lattice to disk
filename=strcat(prefix,'T_',num2str(T(k)),'-.rtable');
dlmwrite(filename,l,' ');
end

% write Cv, S, T to disk
f=[T',Cv',s'];
filename=strcat(prefix,'-TCvs','.rtable');
dlmwrite(filename,f,' ');