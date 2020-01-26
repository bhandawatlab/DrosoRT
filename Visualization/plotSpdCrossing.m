function [] = plotSpdCrossing(condition,opts)
% This function plots the speed when crossing the light border.
%
% Inputs:
%    condition: 'Out2In' if entering the light region; 'In2Out' if leaving
%       the light region
%    opts: structure
%       opts.dataConsFold: folder with the consolidated data
%       opts.crossFold: folder with the crossing tracks data
%       opts.plotFig: whether or not to plot the figure

close all

% set up datafolders and to look in
DataFold = opts.dataConsFold;
crossFold = opts.crossFold;
datainfo = dir(DataFold);
datainfo2 = dir(crossFold);
dataFiles = cell(1,length(datainfo)-2);
dataCrossFiles = cell(1,length(datainfo)-2);

% determine the condition specified by user
if strcmpi(condition,'Out2In')
    sce = 'speed entering';
else
    sce = 'speed leaving';
end

% extract data file names and locations
genotype = struct('name',[],'files',[]);
for K = 1:length(dataFiles)
    thisdir = datainfo(K+2).name;
    thisdir2 = datainfo2((K.*2)+1).name;
    dataFiles{K} = [DataFold '\' thisdir];
    dataCrossFiles{K} = [crossFold '\' thisdir2];
    genotype.name{K} = thisdir(1:end-10);
end

% initialize cel arrays
nGen = length(dataFiles);
allFlies = cell(size(dataFiles));
allTmpCond = cell(size(dataFiles));
allTmpTracksCond = cell(size(dataFiles));
% loop through each genotype in the folder
for i = 1:nGen
    % load the fly data structures for fly tracks
    load(dataFiles{i},'Data','empFlys','curvPks','stopCond','fs')    
    % load the crossing track matrices based on if crossin in or out
    if strcmpi(condition,'Out2In')
        load(dataCrossFiles{i},'tmpOutIn','tmpTracksOutIn2');
        allTmpCond{i} = tmpOutIn;
        allTmpTracksCond{i} = tmpTracksOutIn2;
    else
        load(dataCrossFiles{i},'tmpInOut','tmpTracksInOut2');
        allTmpCond{i} = tmpInOut;
        allTmpTracksCond{i} = tmpTracksInOut2;
    end
    
    % if there are flies that do not cross, then remove data from data
    % structure based on euclidean distance to flies that do cross.
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
end

allSpd = cell(1,nGen);
c = 'gr';m1 = 0;m2 = 0;
% plot tracks from 15 frames before crossing to 30 frames after crossing
nBefore = 15;nAfter = 30;nTot = nBefore+nAfter+1;
% for i = 1:nGen % use this line if interested in all genotypes 
%loop through the first 2 genotypes (control and retinal in this case)
for i = 1:2
    flyN = allTmpCond{i}(:,3)+1;
    first2turns.Bef{i} = [];first2turns.Aft{i} = [];
    laterTurns.Bef{i} = [];laterTurns.Aft{i} = [];
    for j = 1:size(allTmpCond{i},1)
        tmpTrack = allTmpTracksCond{i}(j,~isnan(allTmpTracksCond{i}(j,:)));
        tmpTrack2 = [tmpTrack(1)-nBefore:tmpTrack(1)-1 tmpTrack];
        tmpTrack = tmpTrack2;
        rTmp = sqrt(allFlies{i}.x(flyN(j),tmpTrack).^2+allFlies{i}.y(flyN(j),tmpTrack).^2);
        spdTmp = sqrt(allFlies{i}.thrust(flyN(j),tmpTrack).^2+allFlies{i}.slip(flyN(j),tmpTrack).^2);
        allSpd{i}(j,1:1200) = nan;
        %if ~any(rTmp>opts.boundary)
            allSpd{i}(j,1:min(length(spdTmp),1200)) = spdTmp(1:min(length(spdTmp),1200));
        %end
    end
    allSpd{i}(sum(allSpd{i},2)==0,:) = [];
    
    % get the limits for the plots rounded up to the nearest multiple of 5
    m1 = ceil(max(m1,max(nanmean(allSpd{i}(:,1:nTot),1)))/5)*5;
    m2 = max(m2,round(max(m1+nanstd(allSpd{i}(:,1:nTot)))/5)*5);
    
    % plot the mean +/- standard deviation
    figure(1);hold on
    shadedErrorBar((-nBefore:nAfter)./fs,nanmean(allSpd{i}(:,1:nTot),1),nanstd(allSpd{i}(:,1:nTot)),'lineprops',c(i));
    xlabel('time (s)');ylabel('mm/s');ylim([0 m2])
    suptitle(sce)
    % plot the mean +/- standard estimate of mean
    figure(2);hold on
    shadedErrorBar((-nBefore:nAfter)./fs,nanmean(allSpd{i}(:,1:nTot),1),nanstd(allSpd{i}(:,1:nTot))./sqrt(size(allSpd{i},1)),'lineprops',c(i));
    xlabel('time (s)');ylabel('mm/s');ylim([0 m1]);
    suptitle(sce)
end
% set labels
figure(1);text(0.25,10,'red=mean Retinal +/- STD','Color','r');
text(0.25,13,'green=mean Control +/- STD','Color','g')
figure(2);text(0.25,5,'red=mean Retinal +/- SEM','Color','r');
text(0.25,8,'green=mean Control +/- SEM','Color','g')
% print the figures
if opts.plotFig
    for i = 1:get(gcf,'Number')
        figure(i)
        print('-painters','-dpsc2','PaperFigures.ps','-loose','-append');
    end
end
end
