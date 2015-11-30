% Use Fock representation
% State Encoding: See 'generateFockOperators.m'
% This script will do 4 things:
%     1. Call generateFockOperators();
%     2. Generate a NOON like state of N bosons
%     3. Perform an adiabatic evolution with no particle loss
%     4. Save simulation result to mat-file
%
% pin - input parameters
% outList - required output var list
% pout - return value, a struct that contains all vars that is mentioned in
%        outList.
%
% parameters that need to set in pin are:
%    c0x,c1x,cx0,cx1
function pout=adiabaticEvolution_Fock_reverse_RUN2(pin,outList)

%% set default value for some parameters
N=10;
J1=0;          J2=5;
Ec1=-1;        Ec2=-1;
c0x=0;c1x=0;cx0=0;cx1=0;
dt=1e-5;       avgTol=1e-3;   devTol=1e-3;
hisLen=16;     q=0.99;
spt=0.1;       relaxT=20;     tmax=220;

%% set parameters
varNames=fieldnames(pin);
for i=1:length(varNames)
    eval([varNames{i},'=pin.(''',varNames{i},''');']);
end

%% init operators
generateFockOperators();

%% generate initial state
psiInit=zeros(Dim,1);
psiInit(nn2k(1,N+1))=c0x;
psiInit(nn2k(2,N))=c1x;
psiInit(nn2k(N+1,1))=cx0;
psiInit(nn2k(N,2))=cx1;
% normalize psiInit
psiInit=psiInit/norm(psiInit);


%% adiabatic evolution (adaptive time step)
% dt=1e-1;      % initial time step
% avgTol=1e-3;  % tolerance of avgErr, equals to A in OneNote
%               % this parameter will be changed to max(0,avgTol-devTol/2)
%               % in the following script
% devTol=1e-3;  % tolerance of devErr, equals to B in OneNote
% hisLen=16;       % Length of history for psi. When step size need to change,
%                 % we will search the history for a suitable psi to change
%                 % step size.
% q=0.99;       % common ratio of mean value calculation, see OneNote
% spt=0.1;        % sample interval
% tmax=420;
% relaxT=20;      % H will not change during the last relaxT time


% allocate recorder
nRecord=ceil(tmax/spt+1);
tList=zeros(1,nRecord);
psiList=zeros(Dim,nRecord);
JList=zeros(1,nRecord);
EcList=zeros(1,nRecord);
avgErrList=zeros(1,nRecord);
devErrList=zeros(1,nRecord);

% record the initial state
rCount=1;
tList(rCount)=0;
psiList(:,rCount)=psiInit;
JList(rCount)=J1;
EcList(rCount)=Ec1;

% allocate var space
norErr=zeros(hisLen,1);  % error of norm of psi
devErr=zeros(hisLen,1);  % equals to D(i) in OneNote
avgErr=zeros(hisLen,1);  % equals to M(i) in OneNote
psi=cell(hisLen,1);
tHis=(tmax+1)*ones(hisLen,1); % To record the time stamp of psi{i}.
                              % For a valid psi{i}, tHis(i) records its t,
                              % if psi{i} is invalid, tHis(i)=tmax+1.

% init psi
H1=sparse(a1'*a2+a2'*a1);
H2=sparse((a2'*a2-a1'*a1)^2);
H=sparse(-J1*H1+Ec1/8*H2);
psi{1}=psiInit;
psi{2}=psiInit-1i*H*dt*psiInit;
tHis(1)=0;
tHis(2)=dt;

% init vars
avgTol=max(0,avgTol-devTol/2);     % calc the A' in OneNote
norErr(1)=1-psi{1}'*psi{1};  
norErr(2)=1-psi{2}'*psi{2};
devErr(1)=5/8*devTol;         % init value of D(i), see OneNote
devErr(2)=q*devErr(1)+2*(1-q)*abs(norErr(2)-norErr(1));
avgErr(1)=(1-q)*norErr(1);
avgErr(2)=q*avgErr(1)+(1-q)*norErr(2);
prv=1;cur=2;nxt=3;    % init pointer

% If avgErr of the first step exceed tenth of its tolerance
% devTol/20, then display a warning and divide dt by 2.
while (abs(avgErr(2))>devTol/20)
    warning('Initial dt too large, divide dt by 2 automatically.');
    dt=dt/2;
    % evolve again and re-calc err
    psi{2}=psiInit-1i*H*dt*psiInit;
    tHis(2)=dt;
    norErr(2)=1-psi{2}'*psi{2};
    devErr(2)=q*devErr(1)+2*(1-q)*abs(norErr(2)-norErr(1));
    avgErr(2)=q*avgErr(1)+(1-q)*norErr(2);
end

% main loop of time evolution
t=dt;
while(t<=tmax)
    % Update J, Ec and H
    J = (J1+(J2-J1)*t/(tmax-relaxT))*(t<tmax-relaxT)+J2*(t>=tmax-relaxT);
    Ec= (Ec1+(Ec2-Ec1)*t/(tmax-relaxT))*(t<tmax-relaxT)+Ec2*(t>=tmax-relaxT);
    H=sparse(-J*H1+Ec/8*H2);
    % Try to run one step using current dt
    psi{nxt}=psi{prv}-2i*H*dt*psi{cur};
    tHis(nxt)=t+dt;
    
    % Update err
    norErr(nxt)=1-norm(psi{nxt})^2;
    avgErr(nxt)=q*avgErr(cur)+(1-q)*norErr(nxt);
    devErr(nxt)=q*devErr(cur)+2*(1-q)*abs(norErr(nxt)-norErr(cur));
    
    % Check whether to change dt, see OneNote for criteria
    tmp=avgTol*t/tmax+devTol/2;
    %tmp=devTol/2/(1-avgTol/(avgTol+devTol/2)*(t/tmax)^3);
    if (abs(avgErr(nxt))>tmp || devErr(nxt)>devTol ...
             || (abs(avgErr(nxt))<tmp/4 && devErr(nxt)<devTol/4) )
    %# Change dt
        % Decide how to change dt: 
        % 00: both too small;   01: avgErr too large; 
        % 10: devErr too large  11: both too large
        changeMode=0;
        if abs(avgErr(nxt))>tmp
            changeMode=changeMode+1;
        end
        if (devErr(nxt)>devTol)
            changeMode=changeMode+2;
        end
        
        % Find best change point
        % only points with tHis(i)<=t can be used
        % we are finding the cur point, so we need a prv point, that's
        % why we add criterion tHis(i)>min(tHis)
        bestCur=cur;
        for i=1:hisLen
            if (tHis(i)<t && tHis(i)>min(tHis))
                if ( (changeMode==0 && abs(norErr(i))>abs(norErr(bestCur)) ) ...
                     || (changeMode>0 && abs(norErr(i))<abs(norErr(bestCur))) )
                    bestCur=i;   
                end
            end
        end
        
        % Backtrack to t=tHis(bestCur)
        t=tHis(bestCur);
        cur=bestCur;
        prv=mod(cur-2+hisLen,hisLen)+1;
        nxt=mod(cur,hisLen)+1;
        tHis(tHis>t)=tmax+1;
        
        % Delete invalid records
        while (tList(rCount)>t)
            rCount=rCount-1;
        end
        
        % Change step size evolve again
          % if err too large, dt -> dt/2
        if (changeMode>0)
            display(['dt -> dt/2 @t=',num2str(t),'  dt=',num2str(dt,'%7.1e')]);
            dt=dt/2;
            % General diff format:
            % f(t+dt2) = (1-a^2)*f(t) + a^2*f(t-dt1) + (1+a)*f'(t)*dt2
            % where a=dt2/dt1, second order approximation.
            psi{nxt}=0.75*psi{cur}+0.25*psi{prv}-1.5i*H*psi{cur}*dt;
            tHis(nxt)=t+dt;
          % if err too small, dt -> 1.1*dt
        else
            display(['dt -> 1.1dt @t=',num2str(t)]);
            dt=dt*1.1;
            psi{nxt}=-0.21*psi{cur}+1.21*psi{prv}-2.1i*H*psi{cur}*dt;
            tHis(nxt)=t+dt;
        end
        
        % Reset avgErr and devErr
        avgErr(nxt)=avgErr(cur);
        devErr(nxt)=devErr(cur);
        tmp=avgTol*t/tmax+devTol/2;
        %tmp=devTol/2/(1-avgTol/(avgTol+devTol/2)*(t/tmax)^3);
        gamma=0.2;  % Offset level, see OneNote
        switch changeMode
            case 0
                avgErr(nxt)=avgErr(nxt)+gamma*sign(avgErr(nxt))*tmp*3/4;
                devErr(nxt)=devErr(nxt)+gamma*devTol*3/4;
            case 1
                avgErr(nxt)=avgErr(nxt)-gamma*sign(avgErr(nxt))*tmp*3/4;
            case 2
                devErr(nxt)=devErr(nxt)-gamma*devTol*3/4;
            case 3
                avgErr(nxt)=avgErr(nxt)-gamma*sign(avgErr(nxt))*tmp*3/4;
                devErr(nxt)=devErr(nxt)-gamma*devTol*3/4;
        end
    end
    %# Do not change dt
    % Record state
    if (t-tList(rCount)>=spt || t+dt>tmax)
        rCount=rCount+1;
        tList(rCount)=t;
        psiList(:,rCount)=psi{cur};
        JList(rCount)=J;
        EcList(rCount)=Ec;
        avgErrList(rCount)=avgErr(cur);
        devErrList(rCount)=devErr(cur);
        % plot
%         plotFockState(psi{cur},N,nn2k);
%         set(gca,'ylim',[0 0.6]);
%         title(['t=',num2str(t,'%7.3f'),...
%                '  J=',num2str(J,'%6.3f'),...
%                '  Ec=',num2str(Ec,'%6.3f'),...
%                '  |norErr|=',num2str(abs(norErr(cur)),'%6.1e'),...
%                '  devErr=',num2str(devErr(cur),'%6.1e'),...
%                '  avgErr=',num2str(avgErr(cur),'%6.1e'),...
%                '  dt=',num2str(dt,'%6.1e')]);
%         pause(0.01);
    end
    
    % Loop pointer and t
    prv=mod(prv,hisLen)+1;
    cur=mod(cur,hisLen)+1;
    nxt=mod(nxt,hisLen)+1;
    t=t+dt;
end

%% Data post-processing
finalNorErr=1-psiList(:,rCount)'*psiList(:,rCount);

% trim records
tList=tList(1:rCount);
psiList=psiList(:,1:rCount);
JList=JList(1:rCount);
EcList=EcList(1:rCount);
avgErrList=avgErrList(1:rCount);
devErrList=devErrList(1:rCount);

%% Set return values
for i=1:length(outList)
    pout.(outList{i})=eval(outList{i});
end

end  % end of function
