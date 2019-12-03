function [] = genGenData(opts)
% This function creates a general file for each genotype with information
% about all flies for that genotype
%
% Inputs:
%    opts: Type of points we want to return (3 options)
%       opts.dataFold = folder with the original data files
%       opts.genFold = folder to save the general files to

% get directory names
dataFold = opts.dataFold;
genFolder = opts.genFold;

% get directory information for each genotype
dirinfo = dir(dataFold);
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories

% get the names of each file in each genotype directory
subdir = cell(1,(length(dirinfo)-2));
genotype = struct('name',[],'files',[]);
nGen = length(subdir);
for K = 1:nGen
    thisdir = dirinfo(K+2).name;
    subdir{K} = dir([dataFold '\' thisdir '\*.mat']);
    genotype.name{K} = thisdir;
    
    for KK = 1:length(subdir{K})
        genotype.files{K,KK} = subdir{K}(KK).name;
    end
end

% main loop to process data
fly = cell(1,nGen);curvError = cell(1,nGen);xyError = cell(1,nGen);
for K = 1:nGen
    nFly = length(subdir{K});
    sAll = cell(1,nFly);sArenaAll = cell(1,nFly);
    for KK = 1:nFly
        fName = genotype.files{K,KK};
        oldFolder = [dataFold '\' genotype.name{K}];
        load([oldFolder '\' fName],'s','fs','sArena');
        sAll{KK} = s; sArenaAll{KK} = sArena;
        
        % calculate the curvature
        [curv,xErr,yErr,err] = calcCurvature(s,fs,sArena);
        
        % assign the values to cells
        curvError{K}(KK) = err;
        xyError{K}(:,KK) = [xErr;yErr];
        fly{K}.curvature{1,KK} = curv;
    end
    
    % create file with all the track information for flies within the 
    % current genotype
    createDataGen(sAll,sArenaAll,fly,genotype,curvError,xyError,fs,K,genFolder)
end
end

function [curv,xErr,yErr,err] = calcCurvature(s,fs,sArena)
        
% check if we can create x,y,ang from the thrust, slip, yaw values.
[xErr,yErr] = TSYtoXY(s,fs,sArena);

% create curvature from xy positions
[curv,~] = CalcCurvature(s);

% check if we can create the true curvature using thrust, slip, yaw
[~,err,~] = CalcCurvatureTSY(s,fs,sArena);

end

function [] = createDataGen(sAll,sArenaAll,fly,genotype,curvError,xyError,fs,K,genFolder)
% This just an assignment block
for KK = 1:length(sAll)
    s = sAll{KK};sArena = sArenaAll{KK};
    Data.xHead(KK,1:length(s.Head.x)) = s.Head.x;
    Data.yHead(KK,1:length(s.Head.y)) = s.Head.y;
    Data.x(KK,1:length(s.Center.x)) = s.Center.x;
    Data.y(KK,1:length(s.Center.y)) = s.Center.y;
    Data.thrust(KK,1:length(s.Kinematics.thrust)) = s.Kinematics.thrust;
    Data.slip(KK,1:length(s.Kinematics.slip)) = s.Kinematics.slip;
    Data.yaw(KK,1:length(s.Kinematics.yaw)) = s.Kinematics.yaw;
    Data.curv(KK,1:length(fly{K}.curvature{KK})) = fly{K}.curvature{KK};
    Data.ang(KK,1:length(s.AngVec)) = s.AngVec;
    Data.lightOn(KK) = s.LightOn;

    Arena.center(:,KK) = sArena.arenaCent';
    Arena.rad(KK) = sArena.rad;
    Arena.cF(KK) = sArena.cF;

    gen.files = genotype.files(K,:);
    gen.files = gen.files(~cellfun('isempty',gen.files));
end

gen.name = genotype.name{K};
Error.curv = curvError{K};
Error.xy = xyError{K};

% save the data
save([genFolder '\' genotype.name{K} '_Nov13.mat'],'Data','Arena','fs','Error','gen','-v7.3');

end







