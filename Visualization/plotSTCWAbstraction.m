function [] = plotSTCWAbstraction(opts,GlobMin)
% This function plots sample abstraction on real tracks for sharp
% turn/curved walk delineation model
%
% Inputs:
%    opts: structure with fields
%       opts.dataFold: folder with empirical data by fly
%       opts.stopThresh: speed threshold for stopping
%       opts.boundary: threshold for boundary of the arena
%    GlobMin: best fit parameters for delineating between sharp turns and
%       curved walks

DataFold = opts.dataFold;

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

c = 1;
for K = 1:length(subdir)
    if ~exist(['Figures\' genotype.name{K}], 'dir')
        mkdir(['Figures\' genotype.name{K}])
    end
    
    stopCond = [];boundCond = [];
    
    for KK = 1:length(subdir{K})
        s = sAll{K,KK};
        curv = curvAll{K,KK};
        
        spd = sqrt(s.Kinematics.thrust.^2+s.Kinematics.slip.^2);
        r = sqrt(s.Center.x.^2+s.Center.y.^2);
        stop = spd<opts.stopThresh;
        boundary = r(1:end-1)>opts.boundary;
        stop(boundary) = false;
        
        [startNdx,endNdx,type] = startEndSeq(stop|boundary);
        startNdx = startNdx(type == 0);
        endNdx = endNdx(type == 0);
        shortTracks = find(endNdx-startNdx+1<3);
        if sum(shortTracks>0)
            for KKK = 1:length(shortTracks)
                if stop(startNdx(shortTracks(KKK))-1)
                    stop(startNdx(shortTracks(KKK)):endNdx(shortTracks(KKK)))=true;
                else
                    boundary(startNdx(shortTracks(KKK)):endNdx(shortTracks(KKK)))=true;
                end
            end
        end
        startNdx(shortTracks) = [];
        endNdx(shortTracks) = [];
        
        [stopCond] = computeStopBorder(stop,s,curv,KK,stopCond);
        [boundCond] = computeStopBorder(boundary,s,curv,KK,boundCond);
        for KKK = 1:length(startNdx)
            s2{c}.x = s.Center.x(startNdx(KKK):endNdx(KKK)+1);
            s2{c}.y = s.Center.y(startNdx(KKK):endNdx(KKK)+1);
            curv2{c} = curv(startNdx(KKK):endNdx(KKK));
            c = c+1;
        end
    end
end

tracksToCons = floor(rand(80,1).*length(curv2))+1;
for m = 1:length(tracksToCons)
    c = tracksToCons(m);
    
    if mod(m,20)==1
        figure;k = 1;set(gcf,'Position',[2 42 638 600])
    end
    
    pltInf = [5,4,k];
    
    [Cstarts,Cstops,~] = FindCurvStartStopP(curv2{c},GlobMin);
    shortST = Cstops-Cstarts<2;
    Cstarts(shortST) = [];Cstops(shortST) = [];
    [xModel,yModel,ST,err2] = genAbstractionTrack(s2{c},curv2{c},Cstarts,Cstops);
    
    if ~isempty(pltInf)
        
        subplot(pltInf(1),pltInf(2),pltInf(3));
        plot(s2{c}.x,s2{c}.y,'g','LineWidth',2);hold on
        for i = 1:length(ST.Cstarts)
            plot(s2{c}.x(ST.Cstarts(i):ST.Cstops(i)),s2{c}.y(ST.Cstarts(i):ST.Cstops(i)),'r','LineWidth',2)
        end
        plot(xModel,yModel,'k','LineWidth',1)
        distance = sum(sqrt(diff(s2{c}.x).^2+diff(s2{c}.y).^2));
        errAll = err2./distance.*100;
        title(num2str(err2))
        axis equal
        %axis([xmin xmax ymin ymax])
    end
    k = k+1;
end

end