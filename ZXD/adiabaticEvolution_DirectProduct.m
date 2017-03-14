%Author: George-Gate
%Date: 2016/11/16
%Last Modify Date: 2016/11/17
%--------------------------------------------------------------------------
% Use Direct Product representation
%
% [State Encoding] 
%      Consider N spin-1/2 bosons, there are two DOFs: spin DOF and spatial
%  DOF. So the state of one boson can be expressed as 
%                  psi_single=kron(psi_s, psi_x)                    (1)
%  where psi_s is the 2-dim vector represents the spin state and psi_x is
%  a M-dim vector represents the spatial state of the boson.
%      The state for N bosons is
%
%       psi=kron(psi_1,kron(psi_2,kron(psi_3,... kron(psi_(N-1),psi_N)))...)
%
%  where psi_k is the single particle state in (1) for the k-th boson.
%
%
% This script will do 4 things:
%     1. Generate necessary operators
%     2. Generate a SCS of N bosons
%     3. Perform an adiabatic evolution with no particle loss
%     4. Save simulation result to mat-file
%
%% parameters
N=3;M=4;
J1=1;J2=0;     % J will change from J1 to J2 linearly in time period [0,tmax-relaxT]
Ec1=0;Ec2=-1;  % and remain fixed in [tmax-relaxT, tmax]. The same for Ec.

Dim=(2*M)^N;   % The dimension of the state of N bosons

%% constant operators
sigma_x=[0,  1; 1, 0];    % pauli matrices
sigma_y=[0,-1i;1i, 0];
sigma_z=[1,  0; 0,-1];
I_s0=[1,0;0,1];           % Identity operator for single particle spin DOF
I_x0=diag(ones(M,1));     % Identity operator for single particle spatial DOF
sx0=kron(sigma_x/2,I_x0); % spin operators for single spin-1/2, do nothing on the spatial DOF
sz0=kron(sigma_z/2,I_x0);
I0=kron(I_s0,I_x0);       % Identity operator for single particle
% generate culmulative spin operator Sx=sx1+sx2+sx3+...+sxN
IN=1;
Sx=sx0;
Sz=sz0;
for i=2:N
    IN=sparse(kron(I0,IN));
    Sx=sparse(kron(Sx,I0)+kron(IN,sx0));
    Sz=sparse(kron(Sz,I0)+kron(IN,sz0));
end
clear IN;
Sz_2=Sz*Sz;

%% generate SCS as initial state
psi_s0=[1;sign(J1)]/sqrt(2);     % spin state for single boson, J1 should not be zero
psi_x0=zeros(M,1);psi_x0(1)=1;   % spatial state for single boson
psi0=kron(psi_s0,psi_x0);        % state for single boson
psiSCS=psi0;
for i=2:N
    psiSCS=sparse(kron(psi0,psiSCS));
end

[~,SCSbase]=plotSpinState(psiSCS,N,M,[]);
pause(0.01);

%displayFockState(k2nn,psiSCS,'|psiSCS> = ');
%plotFockState(psiSCS,N,nn2k);


%% adiabatic evolution (adaptive time step)
dt=1e-1;      % initial time step
avgTol=1e-5;  % tolerance of avgErr, equals to A in OneNote
              % this parameter will be changed to max(0,avgTol-devTol/2)
              % in the following script
devTol=1e-5;  % tolerance of devErr, equals to B in OneNote
hisLen=16;       % Length of history for psi. When step size need to change,
                % we will search the history for a suitable psi to change
                % step size.
q=0.99;       % common ratio of mean value calculation, see OneNote
spt=0.1;        % sample interval
tmax=60;
relaxT=20;      % H will not change during the last relaxT time


% parameter sweep
% for tmax=[100:20:260,270:10:400,500:100:1000]
tic;

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
psiList(:,rCount)=psiSCS;
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
H=sparse(-2*sqrt(2)*J1*Sx+Ec1/2*Sz_2);
psi{1}=psiSCS;
psi{2}=psiSCS-1i*H*dt*psiSCS;
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
    psi{2}=psiSCS-1i*H*dt*psiSCS;
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
    H=sparse(-2*sqrt(2)*J*Sx+Ec/2*Sz_2);
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
        [~,SCSbase]=plotSpinState(psi{cur},N,M,SCSbase,gca);
        title(['t=',num2str(t,'%7.3f'),...
               '  J=',num2str(J,'%6.3f'),...
               '  Ec=',num2str(Ec,'%6.3f'),...
               '  |norErr|=',num2str(abs(norErr(cur)),'%6.1e'),...
               '  devErr=',num2str(devErr(cur),'%6.1e'),...
               '  avgErr=',num2str(avgErr(cur),'%6.1e'),...
               '  dt=',num2str(dt,'%6.1e')]);
        pause(0.01);
        % end of plot
    end
    
    % Loop pointer and t
    prv=mod(prv,hisLen)+1;
    cur=mod(cur,hisLen)+1;
    nxt=mod(nxt,hisLen)+1;
    t=t+dt;
end

%% Data post-processing & save to file
finalNorErr=1-psiList(:,rCount)'*psiList(:,rCount);

% trim records
tList=tList(1:rCount);
psiList=psiList(:,1:rCount);
JList=JList(1:rCount);
EcList=EcList(1:rCount);
avgErrList=avgErrList(1:rCount);
devErrList=devErrList(1:rCount);

% save to file
if (~exist('mats\dProduct','dir'))
    mkdir('mats\dProduct');
end
save(['mats\dProduct\tmax=',num2str(tmax),' N=',num2str(N),...
      ' M=',num2str(M),'.mat'],...
    'N','Dim','M','Sx','Sz',...
    'avgTol','devTol','hisLen','q',...
    'tmax','spt','relaxT',...
    'rCount','tList','psiList','JList','EcList',...
    'avgErrList','devErrList','finalNorErr');
%%
% final plot 
figure('name',['Result Snapshot(tmax=',num2str(tmax),'  N=',num2str(N),')'],'position',[50 100 1080 400]);
subplot(1,2,1);
[~,SCSbase]=plotSpinState(psiSCS,N,M,SCSbase,gca);
title('Initial State');
subplot(1,2,2);
[~,SCSbase]=plotSpinState(psiList(:,rCount),N,M,SCSbase,gca);
title('Final State');

display(['tmax=',num2str(tmax),'  Finished.']);
% end of parameter sweep
toc;
% end
