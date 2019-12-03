function [fNum] = plottingTracks(flys,totFly,plotStops,title,fName)
% This function plots the trajectories of flies
%
% Inputs:
%    flys: structure with fields
%       Flys.x, Flys.y: x and y positions
%       Flys.stopsBefore, Flys.stopsDuring: locations of "long stops" to
%           visualize
%       Flys.firstEntry: time of first entry after light is turned on
%    totFly: total number of fly tracks
%    plotStops: true/false of plotting "long stops"
%    title: title/genotype of flies
%    fName: file to print to

close all
firstEntry = flys.firstEntry;fNum = 0;
for i = 1:min(totFly,length(firstEntry))
    if mod(i,20)==1
        k = 1;fNum = fNum+1;
        figure(fNum);%set(gcf,'Position',[9 49 784 918])
        set(gcf,'Position',[1.6667, 41.6667, 638.6667, 599.3333])
    end
    subplot(5,4,k);
    plot(flys.x(i,1:firstEntry(i)),flys.y(i,1:firstEntry(i)),'g');hold on
    plot(flys.x(i,firstEntry(i):end),flys.y(i,firstEntry(i):end),'r');
    
    if plotStops
        plot(flys.stopsBefore{i}(1,:),flys.stopsBefore{i}(2,:),'k*')
        plot(flys.stopsDuring{i}(1,:),flys.stopsDuring{i}(2,:),'b*')
    end
	plotCircle([0,0],1/4,50,'c');
    plotCircle([0,0],1.3/4,50,'c');
    plotCircle([0,0],4/4,50,'k');
    
    axis([-1 1 -1 1])
    axis off
    hold on
    k = k+1;
end
for i = 1:fNum
    figure(i)
    suptitle(title)
    if ~isempty(fName)
        print('-painters','-dpsc2',[fName '.ps'],'-loose','-append');
    end
end

end
