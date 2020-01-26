function [] = plotFirstTurnAnalysis(condition,opts)
% This function plots analysis figures for the first two turns
%
% Inputs:
%    condition: crossing 'In2Out' or 'Out2In'
%    opts: structure with fields
%       opts.dataConsFold: folder with the consolidated data
%       opts.crossFold: folder with the crossing tracks data
%       opts.boundary: threshold for boundary of the arena
%       opts.plotFig: true/false
close all

DataFold = opts.dataConsFold;
crossFold = opts.crossFold;
datainfo = dir(DataFold);
datainfo2 = dir(crossFold);
dataFiles = cell(1,length(datainfo)-2);
dataCrossFiles = cell(1,length(datainfo)-2);

if strcmpi(condition,'Out2In')
    sce = 'when entering';
else
    sce = 'when leaving';
end

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
allTmpCond = cell(size(dataFiles));
allTmpTracksCond = cell(size(dataFiles));
for i = 1:nGen
    load(dataFiles{i},'Data','empFlys','curvPks','stopCond')
    load(dataCrossFiles{i},'tmpInOut','tmpTracksInOut2','tmpOutIn','tmpTracksOutIn2');
    if strcmpi(condition,'Out2In')
        load(dataCrossFiles{i},'tmpOutIn','tmpTracksOutIn2');
        allTmpCond{i} = tmpOutIn;
        allTmpTracksCond{i} = tmpTracksOutIn2;
    else
        load(dataCrossFiles{i},'tmpInOut','tmpTracksInOut2');
        allTmpCond{i} = tmpInOut;
        allTmpTracksCond{i} = tmpTracksInOut2;
    end
    
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
    
    
%     allTmpInOut{i} = tmpOutIn;
%     allTmpTracksInOut{i} = tmpTracksOutIn2;
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
turn.tot = nan(nGen,max(cellfun(@length,allTmpCond)),n2Cons);
turn.initDir = nan(nGen,max(cellfun(@length,allTmpCond)),n2Cons,3);
turn.leaveDir = nan(nGen,max(cellfun(@length,allTmpCond)),n2Cons,3);
turn.returnDir = nan(nGen,max(cellfun(@length,allTmpCond)),n2Cons,3);
turn.endDir = nan(nGen,max(cellfun(@length,allTmpCond)),n2Cons,3);
turn.time2Turn = nan(nGen,max(cellfun(@length,allTmpCond)),n2Cons);
turn.dist2Turn = nan(nGen,max(cellfun(@length,allTmpCond)),n2Cons);
turn.spd2Turn = nan(nGen,max(cellfun(@length,allTmpCond)),n2Cons);
turn.radSpd2Turn = nan(nGen,max(cellfun(@length,allTmpCond)),n2Cons);
spd = cell(1,nGen);
for i = 1:nGen
    
    flyN = allTmpCond{i}(:,3)+1;
    first2turns.Bef{i} = [];first2turns.Aft{i} = [];
    laterTurns.Bef{i} = [];laterTurns.Aft{i} = [];
    kk = ones(1,max(flyN));
    for j = 1:size(allTmpCond{i},1)
        tmpTrack = allTmpTracksCond{i}(j,~isnan(allTmpTracksCond{i}(j,:)));
        
        if length(tmpTrack)>30
            spd{i}(j,:) = sqrt(allFlies{i}.thrust(flyN(j),tmpTrack(1)-30:tmpTrack(1)+30).^2+allFlies{i}.slip(flyN(j),tmpTrack(1)-30:tmpTrack(1)+30).^2);
        end
        rTmp = sqrt(allFlies{i}.x(flyN(j),tmpTrack).^2+allFlies{i}.y(flyN(j),tmpTrack).^2);
        if ~any(rTmp>opts.boundary)
            goodTurn = sort(find(tmpTrack(1)<allST{i}.ndx{flyN(j)} & allST{i}.ndx{flyN(j)}<tmpTrack(end)));
            orientRelCent = ang{i}(flyN(j),:).*direct{i}(flyN(j),:);
            pathRelCent = ang2{i}(flyN(j),:).*direct2{i}(flyN(j),:);
            speedAngRelCent = ang3{i}(flyN(j),:).*direct3{i}(flyN(j),:);
            speed = sqrt(diff(allFlies{i}.x(flyN(j),:)).^2+...
                diff(allFlies{i}.y(flyN(j),:)).^2);
            spdRelCent = speed.*cosd(abs(speedAngRelCent));

            tmpTrackType = zeros(size(tmpTrack));

            for k = 1:min(length(goodTurn),n2Cons)
                currSTNdx = allST{i}.all{flyN(j)}{goodTurn(k)}(:,2);
                currSTNdx = currSTNdx(1):min(currSTNdx(end),10799);
                if length(currSTNdx)<3
                    currSTNdx = [currSTNdx(end)-1 currSTNdx currSTNdx(end)+1];
                end

                ndx = currSTNdx-tmpTrack(1)+1;
                ndx(ndx>length(tmpTrackType)) = [];
                ndx(ndx<1) = [];
                tmpTrackType(ndx) = k;

                turn.tot(i,j,k) = allST{i}.tot{flyN(j)}(goodTurn(k));
                turn.leaveDir(i,j,k,:) = orientRelCent(tmpTrack(1:3));
                turn.returnDir(i,j,k,:) = pathRelCent(tmpTrack(end-2:end));
                turn.initDir(i,j,k,:) = orientRelCent(currSTNdx(1:3)-1);
                turn.endDir(i,j,k,:) = orientRelCent(currSTNdx(end-2:end)+1);

                turn.time2Turn(i,j,k) = max(currSTNdx(1)-tmpTrack(1)+1,0)./30;
                turn.dist2Turn(i,j,k) = sum(speed(tmpTrack(1):currSTNdx(1)));
                turn.spd2Turn(i,j,k) = mean(speed(min([tmpTrack(1),currSTNdx(1)]):currSTNdx(1))).*10.*30;% speed in mm/s
                turn.radSpd2Turn(i,j,k) = mean(spdRelCent(min([tmpTrack(1),currSTNdx(1)]):currSTNdx(1))).*10.*30;% speed in mm/s
            end
            AllTracks{i,flyN(j),kk(flyN(j))} = [tmpTrack;tmpTrackType];
            kk(flyN(j)) = kk(flyN(j))+1;
        else
            aaa=1;
        end
    end
end
figNum = 1;

gen = [1:2];

% plot time before first turn
durBefFirstTurn = reshape(turn.time2Turn(gen,:,1)',[],1);
g =  [repmat(genotype.name(1),max(cellfun(@length,allTmpCond)),1);...
    repmat(genotype.name(2),max(cellfun(@length,allTmpCond)),1)];
g(isnan(durBefFirstTurn)) = [];durBefFirstTurn(isnan(durBefFirstTurn)) = [];
g(durBefFirstTurn>10) = [];durBefFirstTurn(durBefFirstTurn>10) = [];
figure;suptitle(['Time Before Entering First Turn (s) ' sce]);
[ss,~,~] = dabest2(durBefFirstTurn,g,'N');
topRow = [num2str(ss.mdCi(1)) '; ' num2str(ss.md) '; ' num2str(ss.mdCi(2))];
title(topRow);

ndx = contains(g,genotype.name(1));
p1 = ranksum(durBefFirstTurn(ndx,1),durBefFirstTurn(~ndx,1));

% plot speed before first turn
spdBefFirstTurn = reshape(turn.spd2Turn(gen,:,1)',[],1);
g =  [repmat(genotype.name(1),max(cellfun(@length,allTmpCond)),1);...
    repmat(genotype.name(2),max(cellfun(@length,allTmpCond)),1)];
g(isnan(spdBefFirstTurn)) = [];spdBefFirstTurn(isnan(spdBefFirstTurn)) = [];
figure;[ss,~,~] = dabest2(spdBefFirstTurn,g,'N');suptitle(['Speed Before Entering First Turn (mm/s) ' sce])
topRow = [num2str(ss.mdCi(1)) '; ' num2str(ss.md) '; ' num2str(ss.mdCi(2))];
title(topRow)

ndx = contains(g,genotype.name(1));
p2 = ranksum(spdBefFirstTurn(ndx,1),spdBefFirstTurn(~ndx,1));

% % plot radial speed before first turn
% radSpdBefFirstTurn = reshape(turn.radSpd2Turn(gen,:,1)',[],1);
% g =  [repmat(genotype.name(1),max(cellfun(@length,allTmpInOut)),1);...
%     repmat(genotype.name(2),max(cellfun(@length,allTmpInOut)),1)];
% g(isnan(radSpdBefFirstTurn)) = [];radSpdBefFirstTurn(isnan(radSpdBefFirstTurn)) = [];
% figure;[ss,~,~] = dabest2(radSpdBefFirstTurn,g,'N');suptitle('Radial speed Before Entering First Turn (mm/s)')
% topRow = [num2str(ss.mdCi(1)) '; ' num2str(ss.md) '; ' num2str(ss.mdCi(2))];
% title(topRow)
% 
% ndx = contains(g,genotype.name(1));
% p3 = ranksum(radSpdBefFirstTurn(ndx,1),radSpdBefFirstTurn(~ndx,1));


figure;set(gcf,'Position',[2 42 838 924])
for j = 1:2
    allorientBefore = squeeze(turn.initDir(j,:,[1,2],2));
    allorientBefore(isnan(allorientBefore)) = [];
    allTurn = squeeze(turn.tot(j,:,[1,2]));
    allTurn(isnan(allTurn)) = [];
    sameDir = sign(allorientBefore)~=sign(allTurn);
    oppDir = sign(allorientBefore)==sign(allTurn);
    
    subplot(2,1,j);
    c = categorical({'Optimal','non-Optimal'});
    bar(c,[sum(sameDir)./length(sameDir),sum(oppDir)./length(oppDir)]);
    ylim([0 1]);ylabel('Probability')
    title([genotype.name{j} ', n=' num2str(numel(oppDir)) ' turns'])
end
suptitle(sce)

if opts.plotFig
    for i = 1:get(gcf,'Number')
        figure(i);
        print('-painters','-dpsc2','PaperFigures.ps','-loose','-append');
    end
end


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



