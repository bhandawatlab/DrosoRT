function [GlobMinX,GlobMinCFit] = getBestFit(opts)
% This function fits 6 parameters to delineate between sharp turns and 
% curved walks
%
% Inputs:
%    opts: Type of points we want to return (3 options)
%       opts.dataFold = folder with the original data files
% 
% Outputs:
%    GlobMinX: The best fit parameters
%    GlobMinCFit: RMSE of the fit

% get data folder
DataFold = opts.dataFold;
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
        [curvTmp,~] = CalcCurvature(s);

        sAll{K,KK} = s;
        curvAll{K,KK} = curvTmp;
        
        close all
    end
end

Ndx = 1;s = {};curv = {};
for K = 1:length(subdir)
    % loop through each file in the genotype folder
    for KK = 1:length(subdir{K})
        sTmp = sAll{K,KK};
        curvTmp = curvAll{K,KK};
        
        % calulate the speed and radial position
        spd = sqrt(sTmp.Kinematics.thrust.^2+sTmp.Kinematics.slip.^2);
        r = sqrt(sTmp.Center.x.^2+sTmp.Center.y.^2);
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
        
        % assign cell arrays of x, y, and curvature trajectories to use in
        % the fitting process
        for KKK = 1:length(startNdx)
            s{Ndx}.x = sTmp.Center.x(startNdx(KKK):endNdx(KKK));
            s{Ndx}.y = sTmp.Center.y(startNdx(KKK):endNdx(KKK));
            curv{Ndx} = curvTmp(startNdx(KKK):endNdx(KKK)-1);
            Ndx = Ndx+1;
        end
        
    end
end
% fit to the trajectories
[GlobMinX,GlobMinCFit] = calcFit(s,curv);
% save the best fit parameters
save('BestFit5','GlobMinX','GlobMinCFit','-v7.3')
end

function [GlobMinX,GlobMinCFit] = calcFit(s,curv)
% do not display plots
pltInfo = [];
% set up objective function
obj = @(x)objFunc(x(1),x(2),x(3),x(4),x(5),x(6),s,curv,pltInfo);
gs = GlobalSearch('Display','iter','NumTrialPoints',2000);

% fit the objective function
tic;x = objFuncWrap(obj,gs);toc

% refit the objective function until either 5 iterations has passed or a
% "low" error has been acheived
j = 0;
while (obj(x)>5 && j<5)
    display(j)
    xTmp = objFuncWrap(obj,gs);
    % keep the new parameters if it fit better than the previous best set
    % of parameters
    if obj(xTmp) < obj(x)
        x = xTmp;
    end
    j = j+1;
end
toc
GlobMinX = x;
% get the fit from the "best" parameter set
[GlobMinCFit,~,~] = obj(x);
end


