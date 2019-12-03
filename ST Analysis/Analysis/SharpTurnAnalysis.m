function [] = SharpTurnAnalysis(features,spdThresh,sbplt)
% This function plots the analysis figures for how flies perform sharp
% turns
% 
% Inputs:
%    yawThresh: threshold for finding peaks in yaw
%    curvThresh: threshold for finding peaks in curvatures
%    fitFile: file with the best fit parameters to delineate between sharp
%       turns and curved walks
% 
% Inputs:
%    features: struction with fields:
%       features.PeakCurvPerc: when does peak in curv occur during ST 
%       features.lenTrack: length of sharp turn
%       features.spdTrack: speed of sharp turn
%       features.curvPeakNum: number of peaks in curvature/ sharp turn
%       features.PeakYawPerc: when does peak in yaw occur during ST 
%       features.yawPeakDelay: when does peak in yaw occur after ST 
%       features.valYaw: yaw during sharp turn
%       features.valCurv: curvature during sharp turn
%       features.valDiffYawCurv: difference in yaw and curvature
%       features.valSlipAng: slip angle during sharp turn
%       features.valDiffSlipAng: difference in slip angle during sharp turn
%       features.percBackwards: percent of sharp turns that move backwards
%       features.distance: distance traveled during sharp turn
%    spdThresh: threshold for speed of sharp turns
%    sbplt: subplot index

PeakCurvPerc = features.PeakCurvPerc;
lenTrack = features.lenTrack;
spdTrack = features.spdTrack;
curvPeakNum = features.curvPeakNum;

PeakYawPerc = features.PeakYawPerc;
yawPeakDelay = features.yawPeakDelay;
valYaw = features.valYaw;
valCurv = features.valCurv;
valDiffYawCurv = features.valDiffYawCurv;
valSlipAng = features.valSlipAng;
valDiffSlipAng = features.valDiffSlipAng;
percBackwards = features.percBackwards;
distance = features.distance;

% remove empty cells
noTracks = lenTrack==0 | lenTrack >30;
[PeakCurvPerc,lenTrack,spdTrack,curvPeakNum,PeakYawPerc,yawPeakDelay,...
    percBackwards,valYaw,valCurv,valSlipAng,valDiffSlipAng,valDiffYawCurv,distance] = removeTracks...
    (PeakCurvPerc,lenTrack,spdTrack,curvPeakNum,PeakYawPerc,yawPeakDelay,...
    percBackwards,valYaw,valCurv,valSlipAng,valDiffSlipAng,valDiffYawCurv,distance,noTracks);
spdTrack = spdTrack.*30.*10;      % speed in mm/s

len = length(lenTrack);

% turns less than 10 frames
frameThresh = 5;
smallTurns = lenTrack<=frameThresh; %| lenTrack>15;

% speed is too slow
slowPeaks = spdTrack<spdThresh(1) | spdTrack>spdThresh(2);

% there are multiple peaks in curvature
multCurvPeaks = curvPeakNum>1;

% the fly is moving backwards
reverseMvmtPeaks = percBackwards>0.2;

% all condisitons combined
badTracks = (smallTurns | slowPeaks | multCurvPeaks | reverseMvmtPeaks);

% remove tracks with all conditions
[PeakCurvPerc,lenTrack,spdTrack,curvPeakNum,PeakYawPerc,yawPeakDelay,...
    percBackwards,valYaw,valCurv,valSlipAng,valDiffSlipAng,valDiffYawCurv,distance] = removeTracks...
    (PeakCurvPerc,lenTrack,spdTrack,curvPeakNum,PeakYawPerc,yawPeakDelay,...
    percBackwards,valYaw,valCurv,valSlipAng,valDiffSlipAng,valDiffYawCurv,distance,badTracks);

% calculate percentage outliers and the yaws with large delays from curv
percOutlier = sum(PeakYawPerc>1 | PeakYawPerc<0)./length(PeakYawPerc).*100;
percLargeDelay = sum(abs(yawPeakDelay)>15)./length(yawPeakDelay).*100;


goodPeakYaw = abs(PeakYawPerc-0.5)<=0.5;

x = 0:0.025:1;
mVal = max(abs(cell2mat(valYaw(goodPeakYaw))));
x2 = [-2:0.5:mVal./2];
thresh = 0;
[~,sYaw] = plotSpec(valYaw(goodPeakYaw),x,x2,thresh,'Yaw (degrees)',false);

x = 0:0.025:1;
% use with slip angle
mVal = max(abs(cell2mat(valSlipAng(goodPeakYaw))));
x2 = linspace(-2,mVal./3.*2,40); % use with slip angle
thresh = 0;
[~,sdiffSlipAng] = plotSpec(valDiffSlipAng(goodPeakYaw),x,x2,thresh,'Slip Angle (degrees)',false);

if ~isempty(sbplt)
    subplot(sbplt(1),sbplt(2),sbplt(3));
else
    figure;set(gcf,'Position',[2 42 638 600])
end
plot((x(2:end)-diff(x(1:2)./2)).*100,abs(sYaw.mean)./(abs(sYaw.mean)+abs(sdiffSlipAng.mean)).*100);hold on;
plot((x(2:end)-diff(x(1:2)./2)).*100,abs(sdiffSlipAng.mean)./(abs(sYaw.mean)+abs(sdiffSlipAng.mean)).*100);
ylim([0 100])
legend({'Mean % contribution of Yaw','Mean % contribution of Slip Angle'},'Location','Best');
xlabel('% of Track');ylabel('% contribution')

if ~isempty(sbplt)
    title(['Speed: ' num2str(spdThresh(1)) ' to ' num2str(spdThresh(2)) ...
        ' mm/s, ' num2str(round(percOutlier)) '% no yaw, ' num2str(length(spdTrack)) ' Turns'])
else
    suptitle(['Speed: ' num2str(spdThresh(1)) ' to ' num2str(spdThresh(2)) ...
        ' mm/s, ' num2str(round(percOutlier)) '% no yaw, ' num2str(length(spdTrack)) ' Turns'])
end


% if plotFig
%     for i = 1:1
%         figure(i);
%         print('-painters','-dpsc2',['Less_than_1mm.ps'],'-loose','-append');
%     end
% end

end


function [PeakCurvPerc,lenTrack,spdTrack,curvPeakNum,PeakYawPerc,yawPeakDelay,...
    percBackwards,valYaw,valCurv,valSlipAng,valDiffSlipAng,valDiffYawCurv,distance] = removeTracks...
    (PeakCurvPerc,lenTrack,spdTrack,curvPeakNum,PeakYawPerc,yawPeakDelay,...
    percBackwards,valYaw,valCurv,valSlipAng,valDiffSlipAng,valDiffYawCurv,distance,pts2Remove)

PeakCurvPerc(pts2Remove) = [];
lenTrack(pts2Remove) = [];
spdTrack(pts2Remove) = [];
curvPeakNum(pts2Remove) = [];
PeakYawPerc(pts2Remove) = [];
yawPeakDelay(pts2Remove) = [];
percBackwards(pts2Remove) = [];
valYaw(pts2Remove) = [];
valCurv(pts2Remove) = [];
valSlipAng(pts2Remove) = [];
valDiffYawCurv(pts2Remove) = [];
valDiffSlipAng(pts2Remove) = [];
distance(pts2Remove) = [];
end


function [y,stats] = plotSpec(valAll,x,x2,thresh,type,plotting)
% reorientate tracks to have peak in value to be positive
for n = 1:length(valAll)
    [~, ndx] = max(abs(valAll{n}(:)));
    sig = sign(valAll{n}(ndx));
    valAll{n} = valAll{n}.*sig;
end

dt = 1./cellfun('length',valAll);
N = cell(1,length(x)-1);
for n = 1:length(valAll)
%    if 1/dt(n)<=20
        for i = 1:length(x)-1
            ndx = ceil(x(i)/dt(n)):ceil(x(i+1)/dt(n));
            ndx = unique(ceil(ndx));
            ndx(ndx==0) = [];
            if ~isempty(ndx)
                N{i} = [N{i} median(valAll{n}(ndx))];
            end
        end
%    end
end

y = zeros(10,length(x2)-1);
for n = 1:length(N)
    tmp = N{n};
    tmp = reshape(tmp(abs(N{n})>=thresh),[],1);
    %pd = fitdist(tmp,'Kernel');
    %y(n,:) = pdf(pd,x2);
    y(n,:) = histcounts(tmp,x2);
    stats.mean(n) = mean(tmp);
    stats.mad(n) = mad(tmp);
    stats.median(n) = median(tmp);
    stats.std(n) = std(tmp);
    stats.sem(n) = std(tmp)./sqrt(length(tmp));
end
y = y./sum(y,2);
if plotting
    figure;subplot(2,2,1);imagesc(y')
    shading interp
    colormap(hot)
    colorbar
    xlabel('% of track');ylabel(type);title('Absolute Val')
    yticklabels = linspace(x2(1),x2(end),5);
    yticks = linspace(1, size(y, 2), numel(yticklabels));
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
    xticklabels = (0:1:size(y,1))./size(y,1).*100+5;
    xticks = (1:1:size(y,1));
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
    
    [X,Y] = meshgrid((x(2:end)-diff(x(1:2)./2)).*100,x2(2:end)-diff(x2(1:2)./2));
    subplot(2,2,2);surf(X,Y,y','FaceColor','interp','EdgeColor','interp');view(110,20)
    xlabel('% of track');ylabel(type);
    set(gcf,'position',[849 100 824 824])
    %colorbar
end
end






















