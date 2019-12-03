function [] = plotSpdCurvEntering(opts)
% This function plots the speed and curvature of flies entering the light
% zone
%
% Inputs:
%    opts: structure with fields
%       opts.dataConsFold: folder with the consolidated data
%       opts.crossFold: folder with the crossing tracks

close all

DataFold = opts.dataConsFold;
crossFold = opts.crossFold;
datainfo = dir(DataFold);
datainfo2 = dir(crossFold);
dataFiles = cell(1,length(datainfo)-2);
dataCrossFiles = cell(1,length(datainfo)-2);

genotype = struct('name',[],'files',[]);
for K = 1:length(dataFiles)
    thisdir = datainfo(K+2).name;
    thisdir2 = datainfo2((K.*2)+1).name;
    dataFiles{K} = [DataFold '\' thisdir];
    dataCrossFiles{K} = [crossFold '\' thisdir2];
    genotype.name{K} = thisdir(1:end-10);
end

nGen = length(dataFiles);
allFlies = cell(size(dataFiles));allST = cell(size(dataFiles));
allTmpInOut = cell(size(dataFiles));
allTmpTracksInOut = cell(size(dataFiles));
for i = 1:nGen
    load(dataFiles{i},'Data','empFlys','curvPks','stopCond')
    load(dataCrossFiles{i},'tmpInOut','tmpTracksInOut2','tmpOutIn','tmpTracksOutIn2');
    
    l1 = length(empFlys.firstEntry);l2 = length(Data.lightOn);
    if l1~=l2
        l3 = min(size(empFlys.x,2),size(Data.x,2));
        d = zeros(l1,l2);
        for j = 1:l1
            for k = 1:l2
                d(j,k) = sqrt(sum((empFlys.x(j,1:l3)-Data.x(k,1:l3)).^2));
            end
        end
        badFly = find(min(d));
        
        Data.xHead(badFly,:) = [];
        Data.yHead(badFly,:) = [];
        Data.x(badFly,:) = [];
        Data.y(badFly,:) = [];
        Data.thrust(badFly,:) = [];
        Data.slip(badFly,:) = [];
        Data.yaw(badFly,:) = [];
        Data.curv(badFly,:) = [];
        Data.ang(badFly,:) = [];
        Data.lightOn(badFly) = [];
    end
    allFlies{i} = Data;
    allST{i} = curvPks;
    for j = 1:length(stopCond.ndx)
        allDec{i}.ndx{j} = [stopCond.ndx{j} curvPks.ndx{j}];
        allDec{i}.tot{j} = [stopCond.tot{j} curvPks.tot{j}];
        allDec{i}.all{j} = [stopCond.all{j} curvPks.all{j}];
        
        [~,N] = sort(allDec{i}.ndx{j});
        allDec{i}.ndx{j} = allDec{i}.ndx{j}(N);
        allDec{i}.tot{j} = allDec{i}.tot{j}(N);
        allDec{i}.all{j} = allDec{i}.all{j}(N);
    end
    
%     allTmpInOut{i} = tmpInOut;
%     allTmpTracksInOut{i} = tmpTracksInOut2;
    allTmpInOut{i} = tmpOutIn;
    allTmpTracksInOut{i} = tmpTracksOutIn2;
end
%allST = allDec;

% calculate angle relative to inside and curvature
ang = cell(1,nGen);direct = cell(1,nGen);curv = cell(1,nGen);
ang2 = cell(1,nGen);direct2 = cell(1,nGen);
ang3 = cell(1,nGen);direct3 = cell(1,nGen);
for i = 1:nGen
    nFlies = size(allFlies{i}.y,1);
    for j = 1:nFlies
        [ang{i}(j,:),direct{i}(j,:)] = computeDirRelCent(allFlies{i}.x(j,:),...
            allFlies{i}.xHead(j,:),allFlies{i}.y(j,:),allFlies{i}.yHead(j,:));
        
        [ang2{i}(j,:),direct2{i}(j,:)] = computeDirRelCent(allFlies{i}.x(j,1:end-2),...
            allFlies{i}.x(j,3:end),allFlies{i}.y(j,1:end-2),allFlies{i}.y(j,3:end));
        
        [ang3{i}(j,:),direct3{i}(j,:)] = computeDirRelCent(allFlies{i}.x(j,1:end-1),...
            allFlies{i}.x(j,2:end),allFlies{i}.y(j,1:end-1),allFlies{i}.y(j,2:end));
    end
    curv{i} = allFlies{i}.curv;
end

n2Cons = 1000;
turn.tot = nan(nGen,max(cellfun(@length,allTmpInOut)),n2Cons);
turn.initDir = nan(nGen,max(cellfun(@length,allTmpInOut)),n2Cons,3);
turn.leaveDir = nan(nGen,max(cellfun(@length,allTmpInOut)),n2Cons,3);
turn.returnDir = nan(nGen,max(cellfun(@length,allTmpInOut)),n2Cons,3);
turn.endDir = nan(nGen,max(cellfun(@length,allTmpInOut)),n2Cons,3);
turn.time2Turn = nan(nGen,max(cellfun(@length,allTmpInOut)),n2Cons);
turn.dist2Turn = nan(nGen,max(cellfun(@length,allTmpInOut)),n2Cons);
turn.spd2Turn = nan(nGen,max(cellfun(@length,allTmpInOut)),n2Cons);
turn.radSpd2Turn = nan(nGen,max(cellfun(@length,allTmpInOut)),n2Cons);
spd = cell(1,nGen);
allSpd = nan(size(allTmpInOut{2},1),1200);
allCurv = nan(size(allTmpInOut{2},1),1200);
for i = 2%1::nGen
    
    flyN = allTmpInOut{i}(:,3)+1;
    first2turns.Bef{i} = [];first2turns.Aft{i} = [];
    laterTurns.Bef{i} = [];laterTurns.Aft{i} = [];
    kk = ones(1,max(flyN));
    for j = 1:size(allTmpInOut{i},1)
        tmpTrack = allTmpTracksInOut{i}(j,~isnan(allTmpTracksInOut{i}(j,:)));
        tmpTrack2 = [tmpTrack(1)-15:tmpTrack(1)-1 tmpTrack];
        tmpTrack = tmpTrack2;
        if length(tmpTrack)>30
            spd{i}(j,:) = sqrt(allFlies{i}.thrust(flyN(j),tmpTrack(1)-30:tmpTrack(1)+30).^2+allFlies{i}.slip(flyN(j),tmpTrack(1)-30:tmpTrack(1)+30).^2);
        end
        rTmp = sqrt(allFlies{i}.x(flyN(j),tmpTrack).^2+allFlies{i}.y(flyN(j),tmpTrack).^2);
        spdTmp = sqrt(allFlies{i}.thrust(flyN(j),tmpTrack).^2+allFlies{i}.slip(flyN(j),tmpTrack).^2);
        
        curvTmp = curv{i}(flyN(j),tmpTrack);
        curvTmp(spdTmp<0.5) = nan;
        
        if ~any(rTmp>opts.boundary)
            
            allSpd(j,1:min(length(spdTmp),1200)) = spdTmp(1:min(length(spdTmp),1200));
            allCurv(j,1:min(length(curvTmp),1200)) = curvTmp(1:min(length(curvTmp),1200));
            
            
        end
    end
    allCurv = allCurv.*180./pi;
    
    figure;subplot(2,1,1);shadedErrorBar((-15:45)./30,nanmean(allSpd(:,1:61),1),nanstd(allSpd(:,1:61)),'lineprops','k');
    xlabel('time (s)');ylabel('mm/s');
    subplot(2,1,2);shadedErrorBar((-15:45)./30,nanmean(allSpd(:,1:61),1),nanstd(allSpd(:,1:61))./sqrt(size(allSpd,1)),'lineprops','k');
    xlabel('time (s)');ylabel('mm/s');suptitle('speed entering')
    
    figure;subplot(2,1,1);shadedErrorBar((-15:45)./30,nanmean(abs(allCurv(:,1:61)),1),nanstd(abs(allCurv(:,1:61))),'lineprops','k');
    xlabel('time (s)');ylabel('degrees');
    subplot(2,1,2);shadedErrorBar((-15:45)./30,nanmean(abs(allCurv(:,1:61)),1),nanstd(abs(allCurv(:,1:61)))./sqrt(size(allCurv,1)),'lineprops','k');
    xlabel('time (s)');ylabel('degrees');suptitle('curv entering')
    
    
    %     shadedErrorBar(x,mean(y,1),std(y),'lineprops','g');
end
print('-painters','-dpsc2','analysisFiguresEntering.ps','-loose','-append');
end


function [ang,dir] = computeDirRelCent(startX,endX,startY,endY)

% compute if the fly is facing more or less towards the center of the arena
u = [startX-endX;startY-endY;zeros(1,length(startY))];
v = [-startX;-startY;zeros(1,length(startY))];

ang = zeros(1,size(u,2));
% calculate the angles between the center facing vector and the
% direction the fly is moving in before the turn and after the turn
for j = 1:size(u,2)
    ang(j) = atan2(norm(cross(u(:,j),v(:,j))),dot(u(:,j),v(:,j))).*180./pi;
end
% calculate the direction (ccw or cw)
w = cross(u,v);dir = sign(w(3,:));
ang = 180-ang;

% figure(1);i = 1;
% plot([startX(i) endX(i)],[startY(i) endY(i)],'k');hold on
% plot([0 startX(i)],[0 startY(i)],'r');hold off
% title([num2str(ang(i).*dir(i))])
% axis equal

end

function [stat] = computeSlidingStats(xVal,yVal,window,noverlap,xrange)

nSlides = ceil(diff(xrange)./(window-noverlap));
for w = 1:nSlides
    currWind = [(w-1)*(window-noverlap)+xrange(1),(w)*window-(w-1)*noverlap+xrange(1)];
    ndx = xVal<currWind(2) & xVal>=currWind(1);
    stat.med(w) = median(yVal(ndx));
    stat.q(w,:) = quantile(yVal(ndx),[0.25 0.75]);
    
    stat.mean(w) = mean(yVal(ndx));
    stat.std(w) = std(yVal(ndx));
    
    stat.xRange(w,:) = currWind;
    stat.x(w) = mean(currWind);
    
    
    stat.n(w,:) = histcounts(yVal(ndx),[-180:20:180])./length(yVal(ndx));
end

p = polyfit(stat.x(~isnan(stat.med)),stat.med(~isnan(stat.med)),1);
stat.medTL = p;

p = polyfit(stat.x(~isnan(stat.mean)),stat.mean(~isnan(stat.mean)),1);
stat.meanTL = p;

end