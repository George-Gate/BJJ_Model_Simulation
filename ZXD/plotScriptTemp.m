%% plot norErr vs.time
% Use together with fix step size evolution program
% Requirement:
%    psiList, rCount

% calc norErr, devErr, avgErr
q=0.99;
norErrTmp=zeros(rCount,1);
devErr=zeros(rCount,1);
avgErr=zeros(rCount,1);
norErrTmp(1)=real(1-psiList(:,1)'*psiList(:,1));
for i=2:rCount
    norErrTmp(i)=real(1-psiList(:,i)'*psiList(:,i));
    avgErr(i)=q*avgErr(i-1)+(1-q)*norErrTmp(i);
    devErr(i)=q*devErr(i-1)+2*(1-q)*abs(norErrTmp(i)-norErrTmp(i-1));
end

% plot
width=100;   % plot range
center=12e4;
% width=(rCount-1)/2;
% center=(rCount+1)/2;
plot(center-width:center+width-1,norErrTmp(center-width:center+width-1),'-o');hold on;
plot(center-width:center+width-1,devErr(center-width:center+width-1),'g');
plot(center-width:center+width-1,avgErr(center-width:center+width-1),'r');hold off;
xlabel('step number');
set(gca,'xlim',[center-width,center+width-1]);
% clear norErrTmp devErr avgErr;

%% compare psi error under different dt
% Requirement: 
%    1 standard file
%    4 sample file
%    dt, psiList, rCount in each mat-file
%    N, tmax, nn2k in standard file
%    N, tmax of all files should be the same
standardFile='mats\noLoss\N=30\(test2)tmax=100  dt=1e-05.mat';
fileList={'mats\noLoss\N=30\(test2)tmax=100  dt=0.002.mat', ...
          'mats\noLoss\N=30\(test2)tmax=100  dt=0.001.mat', ...
          'mats\noLoss\N=30\(test2)tmax=100  dt=0.0005.mat', ...
          'mats\noLoss\N=30\(test2)tmax=100  dt=0.00025.mat'};
standard=load(standardFile,'N','dt','psiList','tmax','rCount','nn2k');

norErr=1-norm(standard.psiList(:,standard.rCount))^2;
figure('name',['Precision comparision. Standard dt=',num2str(standard.dt),...
       '  tmax=',num2str(standard.tmax),...
       '  norErr=',num2str(norErr,'%9.1e')]);
% read 4 samples and compare
for i=1:4
    ha=subplot(2,2,i);
    sample=load(fileList{i},'dt','psiList','rCount');
    diffErr=plotFockState(abs(sample.psiList(:,sample.rCount))-abs(standard.psiList(:,standard.rCount)),...
                  standard.N,standard.nn2k,ha);
    xlabel(['n_1 (diffErr=',num2str(diffErr,'%9.2e'),')']);
    ylabel('[diff of abs(psi)]^2');
    norErr=1-norm(sample.psiList(:,sample.rCount))^2;
    title(['dt=',num2str(sample.dt),'  norErr=',num2str(norErr,'%9.1e')]);
end
clear standardFile fileList standard sample ha norErr;

%% plot state distribution vs. time
% Requirement: psiList, rCount, Nindex
set(pcolor(abs(psiList([Nindex;1],1:rCount))),'EdgeAlpha',0);