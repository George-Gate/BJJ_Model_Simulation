% Starter for RUN1
% Will do following parameter sweep:
% (Here we assume kappa_i=Omega_i)
% [RUN1.1]
% N=25, J1=40, Omega1=0.002, Omega2=0.002
% (This Omega means 0.5 particle loss every 10 unit of time)
%    tmax = [40 60 80 100 120 140 160 180 200 300 400 500]
% 
% [RUN1.2]
% N=25, J1=8, Omega1=0.002, Omega2=0.002
%    tmax = [24 28 32 36 40 44 48 52 56 76 96 116 136 156 176 196 216 416 616]
% (This tmax list is correspond to the tmax list of RUN1.1, for they have
% the same dJ/dt.)
%
% [RUN1.3]
% N=25, J1=40, Omega1=0.002, Omega2=0.001
% (Different decay rate for two modes)
%    tmax = [40 60 80 100 120 140 160 180 200 300 400 500]
%
% [RUN1.4]
% N=25, J1=40, Omega1=0.0002, Omega2=0.0002
% (This Omega means 0.1 particle loss every 10 unit of time)
%    tmax = [40 60 80 100 120 140 160 180 200 300 400 500 600 700 800 900 1000]
% 
% [RUN1.5]
% N=15, J1=25, Omega1=0.002, Omega2=0.002
%    tmax = [32.5 45 57.5 70 82.5 95 107.5 120 132.5 195 257.5 320]
% (This tmax list is correspond to the tmax list of RUN1.1, for they have
% the same dJ/dt.)
%


clear;
%% default parameters
J2=0;     % J will change from J1 to J2 linearly in time period [0,tmax-relaxT]
Ec1=-1;Ec2=-1;  % and remain fixed in [tmax-relaxT, tmax]. The same for Ec.


%% default parameters for evolution
dt0= 1e-5;    % initial step size
relTol=1e-3; % the maximum relative error at t=tmax
maxDt=1e-1;  % maximum step size
spt=0.1;     % sample interval
relaxT=20;

%% generate task list
% task parameters
runName={'RUN1.1','RUN1.2','RUN1.3','RUN1.4','RUN1.5'};
parameterList={struct('N',25,'J1',40,'Omega1',2e-3,'Omega2',2e-3,... % RUN1.1
                 'tmax',[40 60 80 100 120 140 160 180 200 300 400 500]);...
               struct('N',25,'J1',8,'Omega1',2e-3,'Omega2',2e-3,... % RUN1.2
                 'tmax',[24 28 32 36 40 44 48 52 56 76 96 116 136 156 176 196 216 416 616]);...
               struct('N',25,'J1',40,'Omega1',2e-3,'Omega2',1e-3,... % RUN1.3
                 'tmax',[40 60 80 100 120 140 160 180 200 300 400 500]);...
               struct('N',25,'J1',40,'Omega1',2e-4,'Omega2',2e-4,... % RUN1.4
                 'tmax',[40 60 80 100 120 140 160 180 200 300 400 500 600 700 800 900 1000]);...
               struct('N',15,'J1',25,'Omega1',2e-3,'Omega2',2e-3,... % RUN1.5
                 'tmax',[32.5 45 57.5 70 82.5 95 107.5 120 132.5 195 257.5 320])
               };

% generate task list
taskList=cell(100,1);
taskCount=0;

for runID=1:5
    % load parameters
    N=parameterList{runID}.N;
    J1=parameterList{runID}.J1;
    Omega1=parameterList{runID}.Omega1;kappa1=Omega1;
    Omega2=parameterList{runID}.Omega2;kappa2=Omega2;
    tmax=parameterList{runID}.tmax;
    % generate task
    for i=1:length(tmax)
        clear p;
        p.taskName=[runName{runID},' tmax=',num2str(tmax(i))];
        p.N=N;
        p.J1=J1;p.J2=J2;
        p.Ec1=Ec1;p.Ec2=Ec2;
        p.Omega1=Omega1;
        p.Omega2=Omega2;
        p.kappa1=kappa1;
        p.kappa2=kappa2;
        p.tmax=tmax(i);
        p.dt=dt0;
        p.relTol=relTol; 
        p.maxDt=maxDt;
        p.spt=spt;
        p.relaxT=relaxT;
        p.saveFile=['mats\',runName{runID},'\N=',num2str(p.N),'  tmax=',num2str(p.tmax),...
              '  k1=',num2str(p.kappa1,'%6.1e'),'  k2=',num2str(p.kappa2,'%6.1e'),'  relTol=',num2str(p.relTol,'%6.0e'),'.mat'];
        
        taskCount=taskCount+1;
        taskList{taskCount}=p;
    end
end

%% Set var list to record
outList={'N','Dim','nn2k','k2nn',...
         'tmax','spt','relaxT','relTol','maxDt',...
         'Omega1','Omega2','kappa1','kappa2',...
         'rCount','tList','rhoList','JList','EcList','checkPoint',...
         'hCount','tHistory','dtHistory','errHistory'};

%% Run task list
for i=1:taskCount
    p=taskList{i};
    % Check if the path exist
    [pathstr,~,~]=fileparts(p.saveFile);
    if (~exist(pathstr,'dir'))
        mkdir(pathstr);
    end
    
    % Check if the target mat file exist
    % If exist, we assume that this task was finished before and will skip
    % this task.
    if (exist(p.saveFile,'file')~=2)
        display(['Batch task: ',p.taskName]);
        
        batch([
        ...% run
        'results=masterEQEvolutionWithLoss_Fock_RUN1(p,outList);',...
        ...% save to file
        'save(p.saveFile,''results'');']);
    else
        display(['Skip ',p.taskName,' because save file already exist.']);
    end
end
