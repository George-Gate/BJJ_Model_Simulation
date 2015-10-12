% Use Fock representation
% State Encoding: See 'generateFockOperators.m'
% Use master equation to calculate the evolution of system's reduce density
% matrix when particle loss exists.
% See OneNote for master equation. 
% This script will do 4 things:
%     1. Call generateFockOperators()
%     2. Generate a BJJ Ground state as the initial state
%     3. Perform an adiabatic evolution with particle loss using master
%        equation. This is an adaptive step size program.
%     4. Save simulation result to mat-file
%% parameters
N=25;
J1=8;J2=0;     % J will change from J1 to J2 linearly in time period [0,tmax-relaxT]
Ec1=-1;Ec2=-1;  % and remain fixed in [tmax-relaxT, tmax]. The same for Ec.
Omega1=0.05;
Omega2=0.05;
kappa1=0.05;
kappa2=0.05;

%% init operators
generateFockOperators();
%generateProjectionOperators();


%% adiabatic evolution of Master Equation (Adaptive)
dt= 1e-5;    % initial step size
relTol=1e-3; % the maximum relative error at t=tmax
maxDt=1e-2;  % maximum step size
spt=0.1;     % sample interval
tmax=60;
relaxT=20;



% parameter sweep
%for dt=[0.002,0.001,0.0005,0.00025,1e-5]
tic;

% allocate recorder
nRecord=ceil(tmax/spt+1);
tList=zeros(nRecord,1);
rhoList=cell(nRecord,1);
JList=zeros(nRecord,1);
EcList=zeros(nRecord,1);
hCount=0;
tHistory=zeros(tmax*1e4,1);    % record t, dt and err at each time when 
dtHistory=zeros(tmax*1e4,1);   % an err is calculated
errHistory=zeros(tmax*1e4,1);


% calc constant operators
H1=sparse(a1'*a2+a2'*a1);
H2=sparse((a2'*a2-a1'*a1)^2);
On1On2=sparse(Omega1*a1'*a1+Omega2*a2'*a2);
O2=sparse(kappa1*a1'*a1+kappa2*a2'*a2);

% define some functions
calcJ=   @(t)(J1+(J2-J1)*t/(tmax-relaxT))*(t<tmax-relaxT)+J2*(t>=tmax-relaxT);
calcEc=  @(t)(Ec1+(Ec2-Ec1)*t/(tmax-relaxT))*(t<tmax-relaxT)+Ec2*(t>=tmax-relaxT);
calcO1=  @(J,Ec)sparse(-1i*(-J*H1+Ec/8*H2+On1On2));
calcDrho=@(O1,rho)(O1-O2)*rho-rho*(O1+O2)+2*kappa1*a1*rho*a1'+2*kappa2*a2*rho*a2';

% init rho & pointer
H=sparse(-J1*H1+Ec1/8*H2);
rhoSCS=makeState('BJJ Ground','rho',nn2k,N,Dim,H);
O1=calcO1(J1,Ec1);
drho=calcDrho(O1,rhoSCS);
rho={rhoSCS,rhoSCS+dt*drho,zeros(Dim,Dim)};
prv=1;cur=2;nxt=3;

% init check point
% checkPoint records the last trusted time point
% 'Trusted' means that the error of that time point is acceptable
% checkPoint.rho{1:2}=rho{[prv,cur]}
% checkPoint.t is the time stamp of rho{cur}
% checkPooint.dt is the time interval between rho{prv} and rho{cur}
checkPoint.rho{1}=rho{prv};
checkPoint.rho{2}=rho{cur};
checkPoint.dt=dt;
checkPoint.t=dt;

% init recorder
rCount=1;
tList(rCount)=0;
rhoList{rCount}=rho{1};
JList(rCount)=J1;
EcList(rCount)=Ec1;

t=dt;
stepCounter=1;
while(t<=tmax)
    stepCounter=stepCounter-1;
  
    % Update J, Ec, H and rho
    J  = calcJ(t);
    Ec = calcEc(t);
    O1 = calcO1(J,Ec);
    drho=calcDrho(O1,rho{cur});
    rho{nxt}=rho{prv}+2*dt*drho;

    % Check if need to change step size
    if (stepCounter<=0)
        rhoTmp=rho{nxt};   % record rho{nxt} calc by dt @t+dt
        
        % Try evolve use dt2=dt/2
        % General diff format:
        % f(t+dt2) = (1-a^2)*f(t) + a^2*f(t-dt1) + (1+a)*f'(t)*dt2
        % where a=dt2/dt1, second order approximation.
        rho{nxt}=0.75*rho{cur}+0.25*rho{prv}+0.75*dt*drho;
        J  = calcJ(t+dt/2);
        Ec = calcEc(t+dt/2);
        O1 = calcO1(J,Ec);
        drho=calcDrho(O1,rho{nxt});
        rho{nxt}=rho{cur}+dt*drho; % record rho{nxt} calc by dt/2 @t+dt
        
        % Estimate relative error
        err=2*norm(rho{nxt}-rhoTmp,'fro')/norm(rho{nxt},'fro');
        %errRecord=[errRecord,err/dt];dtRecord=[dtRecord,dt];tRecord=[tRecord,t];
        err=err/dt*tmax;
        % Record history
        hCount=hCount+1;
        tHistory(hCount)=t;    % record t, dt and err at each time when 
        dtHistory(hCount)=dt;   % an err is calculated
        errHistory(hCount)=err;  % err is the estimated final relative error
        
        
        % Change step size
        if (err>relTol)       % err exceed tol
            % dt -> dt/2
            dt=dt/2;
            if (dt<eps)
                warning('dt < eps, simulation do not converge. We will save current state and skip this simulation.');
                break;
            end
            % revert state to last check point
            t=checkPoint.t;
            alpha=dt/checkPoint.dt;
            rho{prv}=checkPoint.rho{1};
            rho{cur}=checkPoint.rho{2};
            while (tList(rCount)>t)    % delete records
                rCount=rCount-1;
            end
            % evolve
            J  = calcJ(t);
            Ec = calcEc(t);
            O1 = calcO1(J,Ec);
            drho=calcDrho(O1,rho{cur});
            rho{nxt}=(1-alpha^2)*rho{cur}+alpha^2*rho{prv}+(1+alpha)*dt*drho;
            % Trust new state
            checkPoint.rho{1}=rho{cur};
            checkPoint.rho{2}=rho{nxt};
            checkPoint.dt=dt;
            checkPoint.t=t+dt;
            % set next check point
            stepCounter=1;
        elseif (err<relTol/4 && 1.1*dt<maxDt)   % err too small and dt not reach maxDt
            % Trust current state
            checkPoint.rho{1}=rho{prv};
            checkPoint.rho{2}=rho{cur};
            checkPoint.dt=dt;
            checkPoint.t=t;
            % dt -> 1.1*dt
            dt=dt*1.1;
            % evolve
            alpha=1.1;
            J  = calcJ(t);
            Ec = calcEc(t);
            O1 = calcO1(J,Ec);
            drho=calcDrho(O1,rho{cur});
            rho{nxt}=(1-alpha^2)*rho{cur}+alpha^2*rho{prv}+(1+alpha)*dt*drho;
            % set next check point
            stepCounter=100;
        else
            % Trust current state
            checkPoint.rho{1}=rho{prv};
            checkPoint.rho{2}=rho{cur};
            checkPoint.dt=dt;
            checkPoint.t=t;
            
            stepCounter=200;    % Number of steps between two check point
        end
    end
    
    % record state
    if (t-tList(rCount)>=spt || t+dt>tmax)
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
               '  dt=',num2str(dt,'%6.1e'),...
               '  <N>=',num2str(abs(trace((a1'*a1+a2'*a2)*rho{cur})),'%6.2f'),...
               '  estErr=',num2str(err,'%6.1e')]);
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
tHistory=tHistory(1:hCount);   
dtHistory=dtHistory(1:hCount); 
errHistory=errHistory(1:hCount);

% save to file
save(['mats\withLoss\test\(MEQ)tmax=',num2str(tmax),...
      '  k1=',num2str(kappa1,'%6.1e'),'  k2=',num2str(kappa2,'%6.1e'),'  relTol=',num2str(relTol,'%6.0e'),'.mat'],...
    'N','Dim','nn2k','k2nn',...
    'tmax','spt','relaxT','relTol','maxDt',...
    'Omega1','Omega2','kappa1','kappa2',...
    'rCount','tList','rhoList','JList','EcList','checkPoint',...
    'hCount','tHistory','dtHistory','errHistory');


display(['tmax=',num2str(tmax),'  Finished.']);
% end of parameter sweep
toc;
%end
