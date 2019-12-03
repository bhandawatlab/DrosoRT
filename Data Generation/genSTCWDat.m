function [] = genSTCWDat(opts,GlobMin)
% This function delineate tracks between sharp turns and curved walks
%
% Inputs:
%    opts: Type of points we want to return (3 options)
%       opts.dataFold = folder with the original data files
%       opts.STCWFold = folder to save the sharp turn curved walk data to

% get data folder
DataFold = opts.dataFold;
STCWFold = opts.STCWFold;
dirinfo = dir(DataFold);
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories

% get file names for each genotype folder in data folder
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

% get the data structures (s) and calculate the curvature (curvTmp) for
% each fly and save them to a cell array
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

errAll = cell(1,length(subdir));
for K = 1:length(subdir)
    figNum = 1;%k = 1;
    curvPks = [];curvWalks = [];stopCond = [];boundCond = [];
    err = zeros(length(subdir{K}),1);
    % loop through each file in the genotype folder
    for KK = 1:length(subdir{K})
        s = sAll{K,KK};
        curv = curvAll{K,KK};
        
        % calulate the speed and radial position
        spd = sqrt(s.Kinematics.thrust.^2+s.Kinematics.slip.^2);
        r = sqrt(s.Center.x.^2+s.Center.y.^2);
        % assign stops and boundary conditions
        stop = spd<opts.stopThresh;
        boundary = r(1:end-1)>opts.boundary;
        stop(boundary) = false;
        
        % assign short (<3 frame) tracks not assigned to stop/boundary to
        % stop/boundary
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
        
        % get the stop and boundary condition structures
        [stopCond] = computeStopBorder(stop,s,curv,KK,stopCond);
        [boundCond] = computeStopBorder(boundary,s,curv,KK,boundCond);
        
        % initiate cells
        curvPks.max{KK} = []; curvPks.tot{KK} = []; curvPks.dur{KK} = []; curvPks.avg{KK} = [];
        curvPks.ndx{KK} = []; curvPks.all{KK} = []; curvPks.dirRelCenterBef{KK} = [];
        curvPks.dirRelCenterAft{KK} = [];
        curvWalks.max{KK} = []; curvWalks.tot{KK} = []; curvWalks.dur{KK} = []; curvWalks.avg{KK} = [];
        curvWalks.ndx{KK} = []; curvWalks.all{KK} = []; curvWalks.dirRelCenterBef{KK} = [];
        curvWalks.dirRelCenterAft{KK} = [];
        
        % update the figure conditions
        if mod(KK,10) == 1
            figure(figNum);set(gcf,'Position',[2 42 638 600])
            k = 1;figNum = figNum+1;
        end
        
        % plot the entire fly trajectory
        figure(figNum-1);subplot(5,4,k);
        plot3(s.Center.x,s.Center.y,1:1:length(s.Center.x),'g','LineWidth',0.5);
        hold on;view(2)
        
        % loop through each non boundary/stop track
        for KKK = 1:length(startNdx)
            % initialize x, y, and curvature trajectories to use in the
            % fitting process
            s2{1}.x = s.Center.x(startNdx(KKK):endNdx(KKK)+1);
            s2{1}.y = s.Center.y(startNdx(KKK):endNdx(KKK)+1);
            curv2{1} = curv(startNdx(KKK):endNdx(KKK));
            % give the subplot parameters
            pltInf = [5,4,k];
            % set up objective function
            obj = @(x)objFunc(x(1),x(2),x(3),x(4),x(5),x(6),s2,curv2,pltInf);
            [err(KK,KKK),curvPksTmp,curvWalksTmp] = obj(GlobMin);
            
            % assign fields for sharp turns and curved walks
            for KKKK = 1:length(curvPksTmp.all)
                curvPksTmp.all{KKKK}(:,2) = curvPksTmp.all{KKKK}(:,2)+startNdx(KKK)-1;
            end
            for KKKK = 1:length(curvWalksTmp.all)
                curvWalksTmp.all{KKKK}(:,2) = curvWalksTmp.all{KKKK}(:,2)+startNdx(KKK)-1;
            end
            curvPks.max{KK} = [curvPks.max{KK} curvPksTmp.max];
            curvPks.tot{KK} = [curvPks.tot{KK} curvPksTmp.tot];
            curvPks.dur{KK} = [curvPks.dur{KK} curvPksTmp.dur];
            curvPks.avg{KK} = [curvPks.avg{KK} curvPksTmp.avg];
            curvPks.ndx{KK} = [curvPks.ndx{KK} curvPksTmp.ndx+startNdx(KKK)-1];
            curvPks.all{KK} = [curvPks.all{KK} curvPksTmp.all];
            curvPks.dirRelCenterBef{KK} = [curvPks.dirRelCenterBef{KK},...
                curvPksTmp.dirRelCenterBef];
            curvPks.dirRelCenterAft{KK} = [curvPks.dirRelCenterAft{KK},...
                curvPksTmp.dirRelCenterAft];
            
            % assign to output matrix
            curvWalks.max{KK} = [curvWalks.max{KK} curvWalksTmp.max];
            curvWalks.tot{KK} = [curvWalks.tot{KK} curvWalksTmp.tot];
            curvWalks.dur{KK} = [curvWalks.dur{KK} curvWalksTmp.dur];
            curvWalks.avg{KK} = [curvWalks.avg{KK} curvWalksTmp.avg];
            curvWalks.ndx{KK} = [curvWalks.ndx{KK} curvWalksTmp.ndx+startNdx(KKK)-1];
            curvWalks.all{KK} = [curvWalks.all{KK} curvWalksTmp.all];
            curvWalks.dirRelCenterBef{KK} = [curvWalks.dirRelCenterBef{KK},...
                curvWalksTmp.dirRelCenterBef];
            curvWalks.dirRelCenterAft{KK} = [curvWalks.dirRelCenterAft{KK},...
                curvWalksTmp.dirRelCenterAft];
            
        end
        k = k+2;
    end
    %checkSce(sAll(K,:),curvPks,curvWalks,boundCond,stopCond)
    
    % create an error figure
    figure(figNum);subplot(2,1,1);
    histogram(err);ylabel('Error (% mm/frame)');
    subplot(2,1,2);
    histogram(log10(err));
    xlabel('log Error (% mm/frame)');ylabel('Count')
    suptitle(genotype.name{K})
    set(gcf,'Position',[2 42 638 600])
    
    % print the figures
    if opts.plotFig
        for f = 1:figNum-1
            figure(f);
            suptitle(genotype.name{K})
            print('-painters','-dpsc2','allFigures.ps','-loose','-append');
        end
    end
    
    % save the data
    save([STCWFold '/' genotype.name{K} ' STCW.mat'],'curvPks','curvWalks',...
        'stopCond','boundCond','err');
    close all
    errAll{K} = err;
end
%sum(err,2)./sum(err>0,2);
%figure;histogram(log10(ans))

end