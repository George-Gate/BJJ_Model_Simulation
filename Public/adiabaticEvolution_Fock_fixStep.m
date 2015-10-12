% Use Fock representation
% State Encoding: See 'generateFockOperators.m'
%
%
%% parameters
N=30;
J1=1;J2=0;     % J will change from J1 to J2 linearly in time period [0,tmax-relaxT]
Ec1=0;Ec2=-1;  % and remain fixed in [tmax-relaxT, tmax]. The same for Ec.

%% init operators
generateFockOperators();

%% generate SCS as initial state
psiSCS=zeros(Dim,1);
for i=0:N
    psiSCS(nn2k(i+1,N-i+1))=sign(J1)^i*1/sqrt(factorial(i)*factorial(N-i));
end
psiSCS=psiSCS*sqrt(factorial(N))/sqrt(2)^N;

%displayFockState(k2nn,psiSCS,'|psiSCS> = ');
%plotFockState(psiSCS,N,nn2k);


%% adiabatic evolution (Fix Step)
dt= 0.001;
spt=0.1;     % sample interval
tmax=200;
relaxT=20;


% parameter sweep
%for dt=[0.002,0.001,0.0005,0.00025,1e-5]
tic;

% allocate recorder
nRecord=ceil(tmax/spt+1);
tList=zeros(1,nRecord);
psiList=zeros(Dim,nRecord);
JList=zeros(1,nRecord);
EcList=zeros(1,nRecord);
rCount=1;
tList(rCount)=0;
psiList(:,rCount)=psiSCS;
JList(rCount)=J1;
EcList(rCount)=Ec1;

% init psi & pointer
H1=sparse(a1'*a2+a2'*a1);
H2=sparse((a2'*a2-a1'*a1)^2);
H=sparse(-J1*H1+Ec1/8*H2);
psi={psiSCS,psiSCS-1i*H*dt*psiSCS,zeros(Dim,1)};
prv=1;cur=2;nxt=3;

t=dt;
while(t<=tmax)
    % update J, Ec, H and psi
    J = (J1+(J2-J1)*t/(tmax-relaxT))*(t<tmax-relaxT)+J2*(t>=tmax-relaxT);
    Ec= (Ec1+(Ec2-Ec1)*t/(tmax-relaxT))*(t<tmax-relaxT)+Ec2*(t>=tmax-relaxT);
    H=sparse(-J*H1+Ec/8*H2);
    
    psi{nxt}=psi{prv}-2i*H*dt*psi{cur};

    % record state
    if (t-tList(rCount)>=spt)
        rCount=rCount+1;
        tList(rCount)=t;
        psiList(:,rCount)=psi{cur};
        JList(rCount)=J;
        EcList(rCount)=Ec;
        % plot
%         plotFockState(psi{cur},N,nn2k);
%         set(gca,'ylim',[0 0.6]);
%         title(['t=',num2str(t,'%7.3f'),...
%                '  J=',num2str(J,'%6.3f'),...
%                '  Ec=',num2str(Ec,'%6.3f'),...
%                '  norErr=',num2str(1-psi{cur}'*psi{cur},'%6.3f')]);
%         pause(0.01);
    end
    
    % loop pointer
    prv=cur;cur=nxt;nxt=6-prv-cur;
    t=t+dt;
end

finalNorm=1-psiList(:,rCount)'*psiList(:,rCount);
% save(['mats\noLoss\N=30\(test2)tmax=',num2str(tmax),'  dt=',num2str(dt),'.mat'],...
%     'N','Dim','nn2k','k2nn','dt','tmax','spt','relaxT',...
%     'rCount','tList','psiList','JList','EcList','finalNorm');


display(['tmax=',num2str(tmax),'  Finished.']);
% end of parameter sweep
toc;
%end
