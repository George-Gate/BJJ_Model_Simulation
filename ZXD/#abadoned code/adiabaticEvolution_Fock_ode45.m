% Use Fock representation
% State Encoding: See 'generateFockOperators.m'
%
%
%% parameters
N=50;
J=1;
Ec=0;

%% init operators
generateFockOperators();

%% generate SCS as initial state
psiSCS=zeros(Dim,1);
for i=0:N
    psiSCS(nn2k(i+1,N-i+1))=sign(J)^i*1/sqrt(factorial(i)*factorial(N-i));
end
psiSCS=psiSCS*sqrt(factorial(N))/sqrt(2)^N;

%displayFockState(k2nn,psiSCS,'|psiSCS> = ');
%plotFockState(psiSCS,N,nn2k);


%% adiabatic evolution
spt=0.1;     % sample interval
tmax=80;
relaxT=20;
RelTol=1e-7;
AbsTol=1e-10;

% parameter sweep
% for tmax=[100:20:260,270:10:400,500:100:1000]
% tic;

% init Hailtonian
H1=sparse(a1'*a2+a2'*a1);
H2=sparse((a2'*a2-a1'*a1)^2);

options=odeset('InitialStep',5e-5,'RelTol',RelTol,'AbsTol',AbsTol);
[tList,psiList]=ode45(@(t,y)BJJ_SEQ(t,y,tmax,relaxT,H1,H2),0:spt:tmax,psiSCS,options);
tList=tList.';
psiList=psiList.';
JList = (1-tList/(tmax-relaxT)).*(tList<tmax-relaxT);
EcList=  -tList/(tmax-relaxT).*(tList<tmax-relaxT)-(tList>=tmax-relaxT);
rCount=length(tList);

finalNorm=psiList(:,rCount)'*psiList(:,rCount);
save(['mats\noLoss\N=50\tmax=',num2str(tmax),'(ode45).mat'],...
    'N','Dim','nn2k','k2nn','tmax','spt','relaxT',...
    'rCount','tList','psiList','JList','EcList',...
    'RelTol','AbsTol','finalNorm');


display(['tmax=',num2str(tmax),'  Finished.']);
% end of parameter sweep
% toc;
% end
