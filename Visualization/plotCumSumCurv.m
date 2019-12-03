function [] = plotCumSumCurv(opts)
% This function plots the cumultive sum of curvatures for model predictions
% from the segmentation of empirical tracks to the four kinematic states
%
% Inputs:
%    opts: structure with fields
%       opts.dataFold: folder with the empirical data by fly
%       opts.STCWFold: folder with the sharp turn/curved walk data
%       opts.plotFig: true/false

DataFold = opts.dataFold;
STCWFold = opts.STCWFold;

dirinfo = dir(DataFold);
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories

subdir = cell(1,(length(dirinfo)-2));
genotype = struct('name',[],'files',[]);
for K = 1:length(subdir)
    thisdir = dirinfo(K+2).name;
    subdir{K} = dir([DataFold '\' thisdir '\*.mat']);
    genotype.name{K} = thisdir;
    
    for KK = 1:length(subdir{K})
        genotype.files{K,KK} = subdir{K}(KK).name;
    end
end

sAll = cell(length(subdir),1);curvAll = cell(length(subdir),1);
for K = 1:length(subdir)
    for KK = 1:length(subdir{K})
        fName = genotype.files{K,KK};
        
        % check if we can create thrust, slip, yaw from the x,y,ang values.
        newDatFile = [DataFold '\' genotype.name{K} '\' fName];
        load(newDatFile,'s');
        
        % create curvature from xy positions
        [curv,~] = CalcCurvature(s);

        sAll{K,KK} = s;
        curvAll{K,KK} = curv;
        close all
    end
end

figNum = 0;
err1  = zeros(2,10);err2  = zeros(2,10);err3  = zeros(2,10);err4  = zeros(2,10);
for K = 1:2%length(subdir)
    load([STCWFold '/' genotype.name{K} ' STCW.mat'],'curvPks','curvWalks',...
        'stopCond','boundCond','err');
    
    nFly = length(curvPks.all);
    for KK = 1:nFly
        
        if mod(KK,12)==1
            k = 1;figNum = figNum+1;
            figure(figNum);set(gcf,'Position',[2 42 838 924])
            suptitle(genotype.name{K})
            figure(figNum+100);set(gcf,'Position',[2 42 838 924])
            suptitle([genotype.name{K} ' no walk ST'])
            figure(figNum+200);set(gcf,'Position',[2 42 838 924])
            suptitle([genotype.name{K} ' no stop ST'])
            figure(figNum+300);set(gcf,'Position',[2 42 838 924])
            suptitle([genotype.name{K} ' no curved walk'])
        end
        
        % stopped ST, walking ST, and curved walks
        [synCurv1,type1] = genSynthCurv(curvPks,stopCond,curvWalks,boundCond,KK);
        
        % stopped ST and curved walks
        [synCurv2,type2] = genSynthCurv([],stopCond,curvWalks,boundCond,KK);
        
        % walking ST and curved walks
        [synCurv3,type3] = genSynthCurv(curvPks,[],curvWalks,boundCond,KK);
        
        % stopped ST and walking ST
        [synCurv4,type4] = genSynthCurv(curvPks,stopCond,[],boundCond,KK);
        
        err1(K,KK) = sqrt(mean((cumsum(synCurv1(1:end-1))-cumsum(curvAll{K,KK})).^2));
        err2(K,KK) = sqrt(mean((cumsum(synCurv2(1:end-1))-cumsum(curvAll{K,KK})).^2));
        err3(K,KK) = sqrt(mean((cumsum(synCurv3(1:end-1))-cumsum(curvAll{K,KK})).^2));
        err4(K,KK) = sqrt(mean((cumsum(synCurv4(1:end-1))-cumsum(curvAll{K,KK})).^2));
        
        plottingFigures(curvAll{K,KK},synCurv1,type1,k,genotype.name{K},figNum);
        plottingFigures(curvAll{K,KK},synCurv2,type2,k,genotype.name{K},figNum+100);
        plottingFigures(curvAll{K,KK},synCurv3,type3,k,genotype.name{K},figNum+200);
        plottingFigures(curvAll{K,KK},synCurv4,type4,k,genotype.name{K},figNum+300);
        
        k = k+1;
    end
end
%close all

err1 = err1(err1>0);err2 = err2(err2>0);
err3 = err3(err3>0);err4 = err4(err4>0);
p1 = ranksum(err1,err2);p2 = ranksum(err1,err3);p3 = ranksum(err1,err4);

g = reshape(ones(size(err1))*[1:1:4],[],1);
errAll = [err1;err2;err3;err4];
labels = {'modeled walk','no walking sharp turns',...
    'no stopped sharp turns','no curved walks'};
figure;boxplot(errAll,g,'Labels',labels);hold on
text(1.2,50,['p=' num2str(p1)]);text(2.2,50,['p=' num2str(p2)]);
text(3.2,50,['p=' num2str(p3)])
ylabel('RMS Error (radians)');title(['Orco Contro/Ret n=' num2str(length(err1))])
xtickangle(30)

if opts.plotFig
    print('-painters','-dpsc2','figurePanels.ps','-loose','-append');
    for i = 1:figNum
        figure(i);
        print('-painters','-dpsc2','figurePanels.ps','-loose','-append');
    end
    for i = 1:figNum
        figure(i+100);
        print('-painters','-dpsc2','figurePanels.ps','-loose','-append');
    end
    for i = 1:figNum
        figure(i+200);
        print('-painters','-dpsc2','figurePanels.ps','-loose','-append');
    end
    for i = 1:figNum
        figure(i+300);
        print('-painters','-dpsc2','figurePanels.ps','-loose','-append');
    end
    close all
    dabest2(errAll,cellstr(num2str(g)),'N');
    figure(1);print('-painters','-dpsc2','figurePanels.ps','-loose','-append');
    %figure(2);print('-painters','-dpsc2','figurePanels.ps','-loose','-append');
end

end

function [synCurv,type] = genSynthCurv(curvPks,stopCond,curvWalks,boundCond,KK)
synCurv = zeros(1,10800);type = zeros(1,10800);

if ~isempty(curvPks)
    [synCurv,type] = genSTStopCurvSequence(curvPks,KK,synCurv,type,1);
end
if ~isempty(stopCond)
    [synCurv,type] = genSTStopCurvSequence(stopCond,KK,synCurv,type,2);
end
if ~isempty(curvWalks)
    [synCurv,type] = genCWCurvSequence(curvWalks,KK,synCurv,type,3);
end
if ~isempty(boundCond)
    [synCurv,type] = genCWCurvSequence(boundCond,KK,synCurv,type,4);
end

end

function [synCurv,type] = genSTStopCurvSequence(s,n,synCurv,type,t)
% generate a sequence of synthetic curvature where the fly turns sharp
% turns in one frame and has a steady turn at straight tracks
synCurv(s.ndx{n}) = synCurv(s.ndx{n})+s.tot{n};

for i = 1:length(s.all{n})
    type(s.all{n}{i}(:,2)) = t;
end

end

function [synCurv,type] = genCWCurvSequence(s,n,synCurv,type,t)
lenS = length(s.all{n});
for i = 1:lenS
    tmpSeg = s.all{n}{i}(:,2);
    synCurv(tmpSeg) = synCurv(tmpSeg)+s.avg{n}(i);
    type(tmpSeg) = t;
end
end

function [] = plottingFigures(curvAll,synCurv,type,k,genName,figNum)
[startNdx,endNdx,t] = startEndSeq(type);

co = 'rkgb';
figure(figNum);
%subplot(4,1,k)
subplot(6,2,k)
plot((1:1:length(curvAll))/30,cumsum(curvAll),'k');hold on;
plot((1:1:length(synCurv))/30,cumsum(synCurv),'r')

loc = max([ceil(max(cumsum(synCurv))) ceil(max(cumsum(curvAll)))])+5;
loc2 = min([ceil(min(cumsum(synCurv))) ceil(min(cumsum(curvAll)))]);
for c = 1:4
    startT = startNdx(t==c);
    endT = endNdx(t==c);
    len = endT-startT+1;
    for cc = 1:length(startT)
        plot([startT(cc):endT(cc)]./30,loc.*ones(1,len(cc)),'Color',co(c),'Linewidth',3);
    end
end
xticks([0:60:360])
xlim([0 360])
if k ==1
    text(100/30,loc2+5,'r=ST','fontsize',8)
    text(2000/30,loc2+5,'g=CW','fontsize',8)
    text(4000/30,loc2+5,'k=Stop','fontsize',8)
    text(6000/30,loc2+5,'b=Bound','fontsize',8)
    xlabel('time (s)');ylabel('cumsum curv (rad)')
end
end



