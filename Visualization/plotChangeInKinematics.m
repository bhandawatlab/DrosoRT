function [] = plotChangeInKinematics(opts)
% This function plots changes in kinematics, loopiness of turns, and
% curvature differences between first 2 turns and later turns
%
% Inputs:
% 	opts.plotFig: true/false

close all
load('Orco Raw Vals.mat','rawVals')
fs = 30;
plotKinematicChanges(rawVals,fs)
if opts.plotFig
    for i = 1:get(gcf,'Number')
        figure(i)
        print('-painters','-dpsc2','PaperFigures.ps','-loose','-append');
    end
end
close all

% plot the loopyness of turns
load([pwd '/Data Full/DataTurnBias/Orco Retinal_BorderChoiceAndTurnBias.mat'],'turnAmount')
[~,~] = plotTurnSize(turnAmount,fs,[1,2,3]);
load([pwd '/Data Full/DataTurnBias/Single Antenna Orco_BorderChoiceAndTurnBias.mat'],'turnAmount')
[~,~] = plotTurnSize(turnAmount,fs,[4,5,6]);

if opts.plotFig
    for i = 1:get(gcf,'Number')
        figure(i)
        print('-painters','-dpsc2','PaperFigures.ps','-loose','-append');
    end
end
close all

% Leaving compare turn curvature
load([pwd '/Data Full/DataTurnBias/Orco Retinal_BorderChoiceAndTurnBias.mat'],'turnAmount')
firstTurnsRet = abs(turnAmount.first2TurnsInOut(1,:));
laterTurnsRet = abs(turnAmount.laterTurnsInOut(1,:));
scenario = ' Leaving';
compareCurv(rawVals,firstTurnsRet,laterTurnsRet,scenario)

% Entering compare turn curvature
load([pwd '/Data Full/DataTurnBias/Orco Retinal_BorderChoiceAndTurnBias.mat'],'turnAmount')
firstTurnsRet = abs(turnAmount.first2TurnsOutIn(1,:));
laterTurnsRet = abs(turnAmount.laterTurnsOutIn(1,:));
scenario = ' Entering';
compareCurv(rawVals,firstTurnsRet,laterTurnsRet,scenario)

% % for first two turns
% xBins = [0:10:360];
% NCont = histcounts(firstTurnsCont.*180./pi,xBins);NCont = NCont./sum(NCont);
% NRet = histcounts(firstTurnsRet.*180./pi,xBins);NRet = NRet./sum(NRet);
% NSing = histcounts(firstTurnsSing.*180./pi,xBins);NSing = NSing./sum(NSing);
% figure;plot(xBins(1:end-1)+5,NCont,'k');hold on;plot(xBins(1:end-1)+5,NRet,'r')
% plot(xBins(1:end-1)+5,NSing,'g');xlim([0 360])
% xlabel('Total curvature (degrees)');ylabel('probability');
% title('First 2 sharp turns after leaving')

if opts.plotFig
    for i = 1:get(gcf,'Number')
        figure(i)
        print('-painters','-dpsc2','PaperFigures.ps','-loose','-append');
    end
end

end

function [] = plotKinematicChanges(rawVals,fs)

spdData = cell2mat(rawVals.spd(4:6))'.*fs;
numPts = cellfun(@numel,rawVals.spd(4:6));
g = [repmat({'Before'},numPts(1),1); repmat({'During Out'},numPts(2),1); repmat({'During In'},numPts(3),1)];
figure;dabest2(spdData,g,'N');suptitle(['Run Speed'])

% for comparing curvature
STData = cell2mat(rawVals.sharpTurn(4:6))';% sharp turns are defined by total curvature
numPts = cellfun(@numel,rawVals.sharpTurn(4:6));
g = [repmat({'Before'},numPts(1),1); repmat({'During Out'},numPts(2),1); repmat({'During In'},numPts(3),1)];
figure;dabest2(STData,g,'N');suptitle(['Sharp Turn Curvature'])

% for comparing stop duration
stopDurData = cell2mat(rawVals.dur(1:3))'./fs;
numPts = cellfun(@numel,rawVals.dur(1:3));
g = [repmat({'Before'},numPts(1),1); repmat({'During Out'},numPts(2),1); repmat({'During In'},numPts(3),1)];
figure;dabest2(stopDurData,g,'N');suptitle(['Stop Duration'])

% for comparing run duration
runDurData = cell2mat(rawVals.dur(4:6))'./fs;
numPts = cellfun(@numel,rawVals.dur(4:6));
g = [repmat({'Before'},numPts(1),1); repmat({'During Out'},numPts(2),1); repmat({'During In'},numPts(3),1)];
figure;dabest2(runDurData,g,'N');suptitle(['Run Duration'])

% for comparing run curvature
runCurvData = cell2mat(rawVals.curvWalk(4:6));% curve walks are defined by average curvature
numPts = cellfun(@numel,rawVals.curvWalk(4:6));
g = [repmat({'Before'},numPts(1),1); repmat({'During Out'},numPts(2),1); repmat({'During In'},numPts(3),1)];
figure;dabest2(runCurvData,g,'N');suptitle(['Run Curvature'])


end

function [] = compareCurv(rawVals,firstTurnsRet,laterTurnsRet,scenario)

% for comparing curvature
firstTurnsRet(firstTurnsRet.*180./pi>max(rawVals.sharpTurn{5})) = [];
laterTurnsRet(laterTurnsRet.*180./pi>max(rawVals.sharpTurn{5})) = [];
allTurns = [firstTurnsRet laterTurnsRet];
STData = [cell2mat(rawVals.sharpTurn(4))'; allTurns'.*180./pi; firstTurnsRet'.*180./pi; laterTurnsRet'.*180./pi];
numPts = [cellfun(@numel,rawVals.sharpTurn(4)), numel(allTurns), numel(firstTurnsRet), numel(laterTurnsRet)];
g = [repmat({'Before'},numPts(1),1); repmat({'During Out'},numPts(2),1);...
    repmat({['First 2 turns'  scenario]},numPts(3),1); repmat({['later turns'  scenario]},numPts(4),1)];
figure;[ss,avr,moes] = dabest2(STData,g,'N');suptitle(['Sharp Turn Curvature' scenario])


end

function [allTurns,h] = plotTurnSize(turnAmount,fs,figNum)
allTurns = [turnAmount.first2TurnsInOut,turnAmount.laterTurnsInOut,turnAmount.first2TurnsOutIn,turnAmount.laterTurnsOutIn];
allTurnCell{1} = turnAmount.first2TurnsInOut;titles{1} = 'first2 turns in->out';
allTurnCell{2} = turnAmount.laterTurnsInOut;titles{2} = 'later turns in->out';
allTurnCell{3} = turnAmount.first2TurnsOutIn;titles{3} = 'first2 turns out->in';
allTurnCell{4} = turnAmount.laterTurnsOutIn;titles{4} = 'later turns out->in';
allTurnCell{5} = allTurns;titles{5} = 'all turns';

figure(figNum(1));xl = linspace(0,50,50);
for i = 1:5
    subplot(3,2,i);
    histogram(allTurnCell{i}(3,allTurnCell{i}(1,:)<0).*100,xl,'Normalization','probability')
    hold on;histogram(allTurnCell{i}(3,allTurnCell{i}(1,:)>0).*100,xl,'Normalization','probability')
    legend({'left turn','right turn'});xlabel('Area under curve (cm^2)');title(titles{i})
end

figure(figNum(2));set(gcf,'Position',[2 42 798 774]);h = nan(2,5);
for i = 1:5
    d1 = allTurnCell{i}(3,allTurnCell{i}(1,:)<0).*100;
    d2 = allTurnCell{i}(3,allTurnCell{i}(1,:)>0).*100;
    [h(1,i),h(2,i)] = kstest2(d1,d2,'Tail','smaller');
    [fl,xl] = ecdf(d1);
    [fu,xu] = ecdf(d2);
    
    subplot(3,2,i);hold on;
    plot(xl,fl);
    plot(xu,fu);
    text(20,0.5,['p=' num2str(h(2,i))])
    text(20,0.4,['nLeft=' num2str(numel(d1))])
    text(20,0.3,['nRight=' num2str(numel(d2))])
    xlim([0 30])
    legend({'left turn','right turn'});
    xlabel('Area under curve (mm^2)');ylabel('cdf')
    title([titles{i}])
end

yl = linspace(0,6,50);xl = linspace(0,1.5,50);
figure(figNum(3));
for i = 1:5
    subplot(3,2,i);
    histogram(allTurnCell{i}(2,allTurnCell{i}(1,:)<0)./fs,xl,'Normalization','probability')
    hold on;histogram(allTurnCell{i}(2,allTurnCell{i}(1,:)>0)./fs,xl,'Normalization','probability')
    legend({'left turn','right turn'});xlabel('distance traveled');title(titles{i})
end

end



