function [curvPks,curvWalks,stopCond,boundCond,err] = genSTCWSynth(synthFlys,GlobMin,opts)
% This function generates sharp turn, curved walk, stop, and boundary
% condition structures for synthetic modeled flies
%
% Inputs
%    synthFlys: structure with fields: x,y,thrust,firstEntry
%    GlobMin: parameters for delineating tracks into sharp turn/curved walks
%    opts: structure with fields: border and stopping threshold
% 
% Outputs:
%    curvPks: structure detailing sharp turn condition
%    curvWalks: structure detailing curved walk condition
%    stopCond: structure detailing stop condition
%    boundCond: structure detailing boundary condition
%    err: rmse

nFlys = length(synthFlys.firstEntry);
figNum = 1;
curvPks = [];curvWalks = [];stopCond = [];boundCond = [];
err = zeros(nFlys,1);
for K = 1:nFlys
    s.Center.x = synthFlys.x(K,:);s.Center.y = synthFlys.y(K,:);
    [curv] = CalcCurvature(s);
    
    spd = synthFlys.thrust(K,1:length(s.Center.x)-1).*10;
    r = sqrt(s.Center.x.^2+s.Center.y.^2);
    stop = spd<opts.stopThresh;
    boundary = r(1:end-1)>opts.boundary;
    stop(boundary) = false;
    
    [startNdx,endNdx,type] = startEndSeq(stop|boundary);
    startNdx = startNdx(type == 0);
    endNdx = endNdx(type == 0);
    shortTracks = find(endNdx-startNdx+1<3);
    if sum(shortTracks>0)
        for KK = 1:length(shortTracks)
            if stop(startNdx(shortTracks(KK))-1)
                stop(startNdx(shortTracks(KK)):endNdx(shortTracks(KK)))=true;
            else
                boundary(startNdx(shortTracks(KK)):endNdx(shortTracks(KK)))=true;
            end
        end
    end
    startNdx(shortTracks) = [];
    endNdx(shortTracks) = [];
    
    [stopCond] = computeStopBorder(stop,s,curv,K,stopCond);
    [boundCond] = computeStopBorder(boundary,s,curv,K,boundCond);
    
    curvPks.max{K} = []; curvPks.tot{K} = []; curvPks.avg{K} = [];
    curvPks.ndx{K} = []; curvPks.all{K} = []; curvPks.dirRelCenterBef{K} = [];
    curvPks.dirRelCenterAft{K} = [];
    curvWalks.max{K} = []; curvWalks.tot{K} = []; curvWalks.avg{K} = [];
    curvWalks.ndx{K} = []; curvWalks.all{K} = []; curvWalks.dirRelCenterBef{K} = [];
    curvWalks.dirRelCenterAft{K} = [];
    
    if mod(K,10) == 1
        %figure(figNum);set(gcf,'Position',[2 42 638 600])
        k = 1;figNum = figNum+1;
    end
    
%     figure(figNum-1);subplot(5,4,k);
%     plot3(s.Center.x,s.Center.y,1:1:length(s.Center.x),'g','LineWidth',0.5);
%     hold on;view(2)
    
    for KK = 1:length(startNdx)
        s2{1}.x = s.Center.x(startNdx(KK):endNdx(KK)+1);
        s2{1}.y = s.Center.y(startNdx(KK):endNdx(KK)+1);
        curv2{1} = curv(startNdx(KK):endNdx(KK));
        pltInf = [];%[5,4,k];
        obj = @(x)objFunc(x(1),x(2),x(3),x(4),x(5),s2,curv2,pltInf);
        [err(K,KK),curvPksTmp,curvWalksTmp] = obj(GlobMin);
        
        for KKK = 1:length(curvPksTmp.all)
            curvPksTmp.all{KKK}(:,2) = curvPksTmp.all{KKK}(:,2)+startNdx(KK)-1;
        end
        for KKK = 1:length(curvWalksTmp.all)
            curvWalksTmp.all{KKK}(:,2) = curvWalksTmp.all{KKK}(:,2)+startNdx(KK)-1;
        end
        
        curvPks.max{K} = [curvPks.max{K} curvPksTmp.max];
        curvPks.tot{K} = [curvPks.tot{K} curvPksTmp.tot];
        curvPks.avg{K} = [curvPks.avg{K} curvPksTmp.avg];
        curvPks.ndx{K} = [curvPks.ndx{K} curvPksTmp.ndx+startNdx(KK)-1];
        curvPks.all{K} = [curvPks.all{K} curvPksTmp.all];
        curvPks.dirRelCenterBef{K} = [curvPks.dirRelCenterBef{K},...
            curvPksTmp.dirRelCenterBef];
        curvPks.dirRelCenterAft{K} = [curvPks.dirRelCenterAft{K},...
            curvPksTmp.dirRelCenterAft];
        
        curvWalks.max{K} = [curvWalks.max{K} curvWalksTmp.max];
        curvWalks.tot{K} = [curvWalks.tot{K} curvWalksTmp.tot];
        curvWalks.avg{K} = [curvWalks.avg{K} curvWalksTmp.avg];
        curvWalks.ndx{K} = [curvWalks.ndx{K} curvWalksTmp.ndx+startNdx(KK)-1];
        curvWalks.all{K} = [curvWalks.all{K} curvWalksTmp.all];
        curvWalks.dirRelCenterBef{K} = [curvWalks.dirRelCenterBef{K},...
            curvWalksTmp.dirRelCenterBef];
        curvWalks.dirRelCenterAft{K} = [curvWalks.dirRelCenterAft{K},...
            curvWalksTmp.dirRelCenterAft];
    end
    k = k+2;
%     checkSce(sAll(K,:),curvPks,curvWalks,boundCond,stopCond)
%     save([STCWFold '/' genotype.name{K} ' STCW.mat'],'curvPks','curvWalks',...
%         'stopCond','boundCond','err');
%     close all
%     errAll{K} = err;
end

end