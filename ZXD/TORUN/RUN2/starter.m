% Starter for RUN2
%**********************************************
% Default parameters:
% N=10;
% J1=0;          J2=5;
% Ec1=-1;        Ec2=-1;
% c0x=0;c1x=0;cx0=0;cx1=0;
% dt=1e-5;       avgTol=1e-3;   devTol=1e-3;
% hisLen=16;     q=0.99;
% spt=0.1;       relaxT=20;     tmax=220;
%**********************************************
% Will do following parameter sweep:
% [RUN2.1]
% c0x=1, cx0=exp(1i*linspace(0,2*pi,100));
% 
% [RUN2.2]
% c0x=x, cx0=sqrt(1-x^2);
%    x = linspace(0,1,100);
%
% [RUN2.3]
% c0x=x, cx0=sqrt(1-x^2)*exp(1i*pi/7);
%    x = linspace(0,1,100);
%
% [RUN2.4]
% c0x=1, cx0=x, cx1=sqrt(1-x^2);
%    x = linspace(0,1,100);
% 
% [RUN2.5]
% c0x=1/sqrt(2), cx0=0.5, cx1=0.5*exp(1i*theta);
%    theta = linspace(0,2*pi,100);
%
% [RUN2.6]
% c0x=1, cx0=sqrt(1-0.3^2), cx1=0.3*exp(1i*theta);
%    theta = linspace(0,2*pi,100);
%
% [RUN2.7]
% c0x=1, cx0=sqrt(1-0.3^2), cx1=0.3*exp(1i*pi/7);
%    tmax = [linspace(40,120,41), linspace(140,440,11)];
%
% [RUN2.8]
% c0x=1, cx0=exp(1i*pi/5);
%    tmax = [linspace(40,120,41), linspace(140,440,11)];
%
clear;

%% generate task list
% task parameters
runName={'RUN2.1','RUN2.2','RUN2.3','RUN2.4','RUN2.5','RUN2.6','RUN2.7','RUN2.8'};
parameterList={struct('c0x',linspace(1,1,100),...                   % RUN2.1
                      'cx0',exp(1i*linspace(0,2*pi,100)),...
                      'cx1',linspace(0,0,100));...
               struct('c0x',linspace(0,1,100),...                   % RUN2.2
                      'cx0',sqrt(1-linspace(0,1,100).^2),...
                      'cx1',linspace(0,0,100));...
               struct('c0x',linspace(0,1,100),...                   % RUN2.3
                      'cx0',sqrt(1-linspace(0,1,100).^2)*exp(1i*pi/7),...
                      'cx1',linspace(0,0,100));...
               struct('c0x',linspace(1,1,100),...                   % RUN2.4
                      'cx0',linspace(0,1,100),...
                      'cx1',sqrt(1-linspace(0,1,100).^2));...
               struct('c0x',linspace(1,1,100)/sqrt(2),...           % RUN2.5
                      'cx0',linspace(0.5,0.5,100),...
                      'cx1',0.5*exp(1i*linspace(0,2*pi,100)));...
               struct('c0x',linspace(1,1,100),...                   % RUN2.6
                      'cx0',sqrt(1-0.3^2)*linspace(1,1,100),...
                      'cx1',0.3*exp(1i*linspace(0,2*pi,100)));...
               struct('c0x',1,'cx0',sqrt(1-0.3^2),'cx1',0.3*exp(1i*pi/7),... % RUN2.7
                      'tmax',[linspace(40,120,41), linspace(140,440,11)]);...
               struct('c0x',1,'cx0',exp(1i*pi/5),'cx1',0,...                % RUN2.8
                      'tmax',[linspace(40,120,41), linspace(140,440,11)])       
               };

% generate task list
taskList=cell(1000,1);
taskCount=0;

for runID=1:6    % RUN2.1 ~ RUN2.6
    % load parameters
    c0x=parameterList{runID}.c0x;
    cx0=parameterList{runID}.cx0;
    cx1=parameterList{runID}.cx1;
    % generate task
    for i=1:length(c0x)
        clear p;
        p.taskName=[runName{runID},' (',num2str(i),')'];
        p.c0x=c0x(i);
        p.cx0=cx0(i);
        p.cx1=cx1(i);
        p.saveFile=['mats\',runName{runID},'\i=',num2str(i),'.mat'];
        
        taskCount=taskCount+1;
        taskList{taskCount}=p;
    end
end

for runID=7:8    % RUN2.7 ~ RUN2.8
    % load parameters
    c0x=parameterList{runID}.c0x;
    cx0=parameterList{runID}.cx0;
    cx1=parameterList{runID}.cx1;
    tmax=parameterList{runID}.tmax;
    % generate task
    for i=1:length(tmax)
        clear p;
        p.taskName=[runName{runID},' (',num2str(i),')'];
        p.c0x=c0x;
        p.cx0=cx0;
        p.cx1=cx1;
        p.tmax=tmax(i);
        p.saveFile=['mats\',runName{runID},'\i=',num2str(i),'.mat'];

        taskCount=taskCount+1;
        taskList{taskCount}=p;
    end
end

%% Set var list to record
outList={'N','Dim','nn2k','k2nn','c0x','cx0','cx1'...
         'avgTol','devTol','hisLen','q',...
         'tmax','spt','relaxT',...
         'rCount','tList','psiList','JList','EcList',...
         'avgErrList','devErrList','finalNorErr'};

%% Run task list
pGroup=cell(round(taskCount/5),1);
counter=0;
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
        counter=counter+1;
        pGroup{counter}=p;
    else
        display(['Skip ',p.taskName,' because save file already exist.']);
    end
    
    % Batch many tasks at one time 
    if (counter>=round((taskCount-i+1)/5) || i>=taskCount)
        display(['Batch task: ',pGroup{1}.taskName,' to ',pGroup{counter}.taskName]);
        batch('groupUnpacker(pGroup,counter,outList);');
        counter=0;
    end
    
end
