% Use Fock representation
% State Encoding: See 'generateFockOperators.m'
% Use master equation to calculate the evolution of system's reduce density
% matrix when particle loss exists.
% See OneNote for master equation. 
%% parameters
N=15;
J1=4;J2=0;     % J will change from J1 to J2 linearly in time period [0,tmax-relaxT]
Ec1=-1;Ec2=-1;  % and remain fixed in [tmax-relaxT, tmax]. The same for Ec.
Omega1=0.002;
Omega2=0.002;
kappa1=0.002;
kappa2=0.002;

%% init operators
generateFockOperators();
generateProjectionOperators();


%% adiabatic evolution of Master Equation (Fix Step)
dt= 0.001;
spt=0.1;     % sample interval
tmax=40;
relaxT=20;


% parameter sweep
%for dt=[0.002,0.001,0.0005,0.00025,1e-5]
tic;

% allocate recorder
nRecord=ceil(tmax/spt+1);
tList=zeros(1,nRecord);
rhoList=cell(nRecord,1);
JList=zeros(1,nRecord);
EcList=zeros(1,nRecord);


% calc constant operators
H1=sparse(a1'*a2+a2'*a1);
H2=sparse((a2'*a2-a1'*a1)^2);
On1On2=sparse(Omega1*a1'*a1+Omega2*a2'*a2);
O2=sparse(kappa1*a1'*a1+kappa2*a2'*a2);

% init rho & pointer
H=sparse(-J1*H1+Ec1/8*H2);
rhoSCS=makeState('BJJ Ground','rho',nn2k,N,Dim,H);
O1=sparse(-1i*(-J1*H1+Ec1/8*H2+On1On2));
drho=(O1-O2)*rhoSCS-rhoSCS*(O1+O2)...
        +2*kappa1*a1*rhoSCS*a1'+2*kappa2*a2*rhoSCS*a2';
rho={rhoSCS,rhoSCS+dt*drho,zeros(Dim,Dim)};
prv=1;cur=2;nxt=3;

% init recorder
rCount=1;
tList(rCount)=0;
rhoList{rCount}=rho{1};
JList(rCount)=J1;
EcList(rCount)=Ec1;

t=dt;
while(t<=tmax)
    % update J, Ec, H and rho
    J = (J1+(J2-J1)*t/(tmax-relaxT))*(t<tmax-relaxT)+J2*(t>=tmax-relaxT);
    Ec= (Ec1+(Ec2-Ec1)*t/(tmax-relaxT))*(t<tmax-relaxT)+Ec2*(t>=tmax-relaxT);
    O1=sparse(-1i*(-J*H1+Ec/8*H2+On1On2));
    
    drho=(O1-O2)*rho{cur}-rho{cur}*(O1+O2)...
        +2*kappa1*a1*rho{cur}*a1'+2*kappa2*a2*rho{cur}*a2';
    rho{nxt}=rho{prv}+2*dt*drho;

    % record state
    if (t-tList(rCount)>=spt)
        rCount=rCount+1;
        tList(rCount)=t;
        rhoList{rCount}=rho{cur};
        JList(rCount)=J;
        EcList(rCount)=Ec;
        % plot
        plotFockState(rho{cur},[0 N],nn2k);
        set(gca,'zlim',[0 0.5]);
        title(['t=',num2str(t,'%7.3f'),...
               '  J=',num2str(J,'%6.3f'),...
               '  Ec=',num2str(Ec,'%6.3f'),...
               '  norErr=',num2str(abs(1-trace(rho{cur})),'%6.1e'),...
               '  <N>=',num2str(abs(trace((a1'*a1+a2'*a2)*rho{cur})),'%6.2f')]);
        pause(0.01);
    end

    % loop pointer
    prv=cur;cur=nxt;nxt=6-prv-cur;
    t=t+dt;
end

%% Data post-processing & save to file
finalNorErr=1-trace(rhoList{rCount});

% trim records
tList=tList(1:rCount);
%rhoList=rhoList{1:rCount};  % cannot trim a cell array
JList=JList(1:rCount);
EcList=EcList(1:rCount);

% save to file
save(['mats\withLoss\test\(MEQ)tmax=',num2str(tmax),'  dt=',num2str(dt),'.mat'],...
    'N','Dim','nn2k','k2nn','dt','tmax','spt','relaxT','proj','projNOON',...
    'Omega1','Omega2','kappa1','kappa2',...
    'rCount','tList','rhoList','JList','EcList','finalNorErr');


display(['tmax=',num2str(tmax),'  Finished.']);
% end of parameter sweep
toc;
%end
