% Use Fock representation
% State Encoding: See 'generateFockOperators.m'
% Use master equation to calculate the evolution of system's reduce density
% matrix when particle loss do not exists.
% rho'=-i[H,rho]
%% parameters
N=15;
J1=20;J2=0;     % J will change from J1 to J2 linearly in time period [0,tmax-relaxT]
Ec1=-1;Ec2=-1;  % and remain fixed in [tmax-relaxT, tmax]. The same for Ec.

%% init operators
generateFockOperators();
generateProjectionOperators();


%% adiabatic evolution of Master Equation (Fix Step)
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
rhoList=cell(nRecord,1);
JList=zeros(1,nRecord);
EcList=zeros(1,nRecord);


% init rho & pointer
H1=sparse(a1'*a2+a2'*a1);
H2=sparse((a2'*a2-a1'*a1)^2);
H=sparse(-J1*H1+Ec1/8*H2);
rhoSCS=makeState('BJJ Ground','rho',nn2k,N,Dim,H);
rho={rhoSCS,rhoSCS-1i*dt*(H*rhoSCS-rhoSCS*H),zeros(Dim,Dim)};
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
    H=sparse(-J*H1+Ec/8*H2);

    rho{nxt}=rho{prv}-2i*dt*(H*rho{cur}-rho{cur}*H);

    % record state
    if (t-tList(rCount)>=spt)
        rCount=rCount+1;
        tList(rCount)=t;
        rhoList{rCount}=rho{cur};
        JList(rCount)=J;
        EcList(rCount)=Ec;
        % plot
%         plotFockState(rho{cur},N,proj);
%         set(gca,'ylim',[0 0.6]);
%         title(['t=',num2str(t,'%7.3f'),...
%                '  J=',num2str(J,'%6.3f'),...
%                '  Ec=',num2str(Ec,'%6.3f'),...
%                '  norErr=',num2str(1-trace(rho{cur}),'%6.3f')]);
%         pause(0.01);
    end

    % loop pointer
    prv=cur;cur=nxt;nxt=6-prv-cur;
    t=t+dt;
end

finalNorErr=1-trace(rhoList{rCount});
save(['mats\noLoss\test\(MEQ)tmax=',num2str(tmax),'  dt=',num2str(dt),'.mat'],...
    'N','Dim','nn2k','k2nn','dt','tmax','spt','relaxT','proj','projNOON',...
    'rCount','tList','rhoList','JList','EcList','finalNorErr');


display(['tmax=',num2str(tmax),'  Finished.']);
% end of parameter sweep
toc;
%end
