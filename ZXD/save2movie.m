% Save simulation result as a movie
% [Required Variables]
%  rCount, psiList, tList, JList, EcList, avgErrList, devErrList
%  spt
% [Input]
%  filename
%
%
%%
% create a movie writer
filename='demo(verySlow_evolution).mp4';
videoObj=VideoWriter(filename,'MPEG-4');
videoObj.FrameRate=15*max(1,floor(1/spt));
videoObj.open();

for i=1:rCount
    % load pars
    t=tList(i);
    J=JList(i);
    Ec=EcList(i);
    avgErr=avgErrList(i);
    devErr=devErrList(i);
    % plot
    plotFockState(psiList(:,i),N,nn2k);
    set(gca,'ylim',[0 0.6]);
    title(['t=',num2str(t,'%7.1f'),...
           '  J=',num2str(J,'%6.2f'),...
           '  Ec=',num2str(Ec,'%6.2f'),...
           '  |norErr|=',num2str(abs(1-norm(psiList(:,i))^2),'%6.1e'),...
           '  devErr=',num2str(devErr,'%6.1e'),...
           '  avgErr=',num2str(avgErr,'%6.1e')]);
    pause(0.01);
    % write to file
    videoObj.writeVideo(getframe(gcf));
end

videoObj.close();