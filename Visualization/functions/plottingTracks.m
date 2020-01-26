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
border = 0.985;%1./1.02;

close all
firstEntry = flys.firstEntry;fNum = 0;
for i = 1:min(totFly,length(firstEntry))
    if mod(i,20)==1
        k = 1;fNum = fNum+1;
        figure(fNum);%set(gcf,'Position',[9 49 784 918])
        set(gcf,'Position',[1.6667, 41.6667, 638.6667, 599.3333])
    end
    subplot(5,4,k);hold on
%     plot(flys.x(i,1:firstEntry(i)),flys.y(i,1:firstEntry(i)),'g');hold on
%     plot(flys.x(i,firstEntry(i):end),flys.y(i,firstEntry(i):end),'r');
    tmpR = sqrt(flys.x(i,:).^2+flys.y(i,:).^2);
    [startNdx,endNdx,type] = startEndSeq(tmpR>border);
    startNdx(type==0) = [];endNdx(type==0) = [];
    len = endNdx-startNdx;
    
    plotCircNdx = [startNdx(len<2);endNdx(len<2)];
    pltSeq = zeros(1,size(flys.y,2));
    for j = 1:size(plotCircNdx,2)
        pltSeq(plotCircNdx(:,j)) = 1;
    end
    [startNdx,endNdx,type] = startEndSeq(pltSeq);
    
    
    for j = 1:length(startNdx)
        if type(j) == 1
            if startNdx(j)<firstEntry(i)
                %plot(flys.x(i,startNdx(j):endNdx(j)),flys.y(i,startNdx(j):endNdx(j)),'--g');
                drawArc([0;0],[flys.x(i,startNdx(j)); flys.y(i,startNdx(j))],...
                    [flys.x(i,endNdx(j)); flys.y(i,endNdx(j))],'g') ;
            else
                %plot(flys.x(i,startNdx(j):endNdx(j)),flys.y(i,startNdx(j):endNdx(j)),'--r');
                drawArc([0;0],[flys.x(i,startNdx(j)); flys.y(i,startNdx(j))],...
                    [flys.x(i,endNdx(j)); flys.y(i,endNdx(j))],'r') ;
            end
        else
            if startNdx(j)<firstEntry(i) && endNdx(j)>firstEntry(i)
                plot(flys.x(i,startNdx(j):firstEntry(i)),flys.y(i,startNdx(j):firstEntry(i)),'g');hold on
                plot(flys.x(i,firstEntry(i):endNdx(j)),flys.y(i,firstEntry(i):endNdx(j)),'r');
            elseif startNdx(j)>firstEntry(i)
                plot(flys.x(i,startNdx(j):endNdx(j)),flys.y(i,startNdx(j):endNdx(j)),'r');
            else
                plot(flys.x(i,startNdx(j):endNdx(j)),flys.y(i,startNdx(j):endNdx(j)),'g');
            end
        end
    end
    
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

function [x,y] = drawArc(CentPt,StartPt,EndPt,format)  

r = mean([sqrt(sum((StartPt-CentPt).^2)),sqrt(sum((EndPt-CentPt).^2))]);
thStart = myatan(StartPt(1),StartPt(2),'degrees',2);
thEnd = myatan(EndPt(1),EndPt(2),'degrees',2);

if thStart-thEnd>thEnd-thStart+360
    thEnd = thEnd+360;
end
th = linspace(thStart,thEnd,10);

x = r * cosd(th) + CentPt(1);
y = r * sind(th) + CentPt(2);
plot(x, y,format);
end
