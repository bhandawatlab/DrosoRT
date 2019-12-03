function [features] = processFeatures(yawThresh,curvThresh,opts)
% This function generates features of each sharp turn such as yaw, curv,
% etc from empirical data files and thresholds for yaw and curvature
% 
% Inputs:
%    yawThresh: threshold for finding peaks in yaw
%    curvThresh: threshold for finding peaks in curvatures
%    opts: structure with fields
%       opts.dataFold: folder with the empirical data by flies
%       opts.STCWFold: folder with sharp turn and curved walk information
% 
% Outputs:
%    features: struction with fields:
%       features.PeakCurvPerc: when does peak in curv occur during ST 
%       features.lenTrack: length of sharp turn
%       features.spdTrack: speed of sharp turn
%       features.curvPeakNum: number of peaks in curvature/ sharp turn
%       features.PeakYawPerc: when does peak in yaw occur during ST 
%       features.yawPeakDelay: when does peak in yaw occur after ST 
%       features.valYaw: yaw during sharp turn
%       features.valCurv: curvature during sharp turn
%       features.valDiffYawCurv: difference in yaw and curvature
%       features.valSlipAng: slip angle during sharp turn
%       features.valDiffSlipAng: difference in slip angle during sharp turn
%       features.percBackwards: percent of sharp turns that move backwards
%       features.distance: distance traveled during sharp turn

DataFold = opts.dataFold;
STCWFold = opts.STCWFold;

datainfo = dir(DataFold);
datainfo(~[datainfo.isdir]) = [];  %remove non-directories
STCWinfo = dir(STCWFold);

subdirData = cell(1,(length(datainfo)-2));
genotypeSTCW = struct('files',[]);
genotype = struct('name',[],'files',[]);
for K = 1:length(subdirData)%length(subdirData)-2:length(subdirData)-1
    thisdir = datainfo(K+2).name;
    STCWdir = STCWinfo(K+2).name;
    subdirData{K} = dir([DataFold '\' thisdir '\*.mat']);
    genotypeSTCW.files{K} = STCWdir;
    genotype.name{K} = thisdir;
    
    for KK = 1:length(subdirData{K})
        genotype.files{K,KK} = subdirData{K}(KK).name;
    end
end

nGen =length(subdirData);
nFly = 1;

% generate sharp turns and curved walks
YawPks = cell(nGen,nFly);
sAll = cell(nGen,nFly);
curvAll = cell(nGen,nFly);
sArena = cell(nGen,nFly);
ST = cell(nGen,nFly);
for K = 1:nGen%nGen-1:nGen
    STCWFile = [STCWFold '\' genotypeSTCW.files{K}];
    load(STCWFile,'curvPks');
    for KK = 1:length(subdirData{K})
        fName = genotype.files{K,KK};
        newDatFile = [DataFold '\' genotype.name{K} '\' fName];
        tmp = load(newDatFile,'s','sArena');
        [curv,~] = CalcCurvature(tmp.s);
        sAll{K,KK} = tmp.s;
        curvAll{K,KK} = curv;
        sArena{K,KK} = tmp.sArena;
        currYaw = sAll{K,KK}.Kinematics.yaw.*pi./180;
        [~,pks] = findpeaks(abs(currYaw),'MinPeakHeight',yawThresh);
        YawPks{K,KK} = pks;
        
        ST{K,KK}.ndx = curvPks.ndx{KK};
        ST{K,KK}.all = curvPks.all{KK};
    end
end

sAll = sAll(~cellfun(@isempty,sAll));
curvAll = curvAll(~cellfun(@isempty,curvAll));
sArena = sArena(~cellfun(@isempty,sArena));
YawPks = YawPks(~cellfun(@isempty,YawPks));
ST = ST(~cellfun(@isempty,ST));
nCt = length(ST);

PeakCurvPerc = [];lenTrack = [];spdTrack = [];curvPeakNum = [];
PeakYawPerc = [];yawPeakDelay = [];val = [];percBackwards = [];
valDiffYawCurv = [];distance = [];
for i = 1:nCt
    if ~isempty(ST{i})
        %[Data] = convertToUnitCircle(sAll{i},sArena{i});
        %r = sqrt(Data.Center.x.^2+Data.Center.y.^2);
        spd = sqrt(diff(sAll{i}.Center.x).^2+diff(sAll{i}.Center.y).^2);

        % get kinematics
        currYaw = sAll{i}.Kinematics.yaw;
        currSlip = sAll{i}.Kinematics.slip;
        currCurv = curvAll{i};
        mvmtAng = sAll{i}.AngVec(2:end)-myatan(diff(sAll{i}.Center.x),diff(sAll{i}.Center.y),'degrees',2);
        mvmtAng(mvmtAng>180) = mvmtAng(mvmtAng>180)-floor(mvmtAng(mvmtAng>180)./180).*180;
        mvmtAng(mvmtAng<-180) = mvmtAng(mvmtAng<-180)-ceil(mvmtAng(mvmtAng<-180)./180).*180;
        dmvmtAng = [0 diff(mvmtAng)];

        % obtain delays between yaw peak and the sharp turn
        currYawPks = YawPks{i};
        currNdx = ST{i}.ndx;
        [closestPtBef,closestPtAft,~] = findBeforeAfter(currNdx,currYawPks,'both');
        befDiff = closestPtBef-currNdx;
        aftDiff = closestPtAft-currNdx;
        [~,I] = min([abs(aftDiff);abs(befDiff)],[],1);
        tmp = aftDiff;
        tmp(I==2) = befDiff(I==2);
        [~,ia,~] = intersect(currNdx,currYawPks);
        tmp(ia) = 0;
        yawPeakNdx = currNdx+tmp;


        for j = 1:length(ST{i}.ndx)
            currNdx = ST{i}.ndx(j);
            currTrack = ST{i}.all{j}(:,2);

            if length(currTrack)>2 %&& r(currNdx)<=3.5
                curvTmp = abs(curvAll{i}(currTrack));
                [~,curvPeak] = findpeaks(curvTmp,'MinPeakHeight',curvThresh);

                PeakCurvPerc(i,j) = (currNdx-currTrack(1)+1)./length(currTrack);
                lenTrack(i,j) = length(currTrack);
                spdTrack(i,j) = mean(spd(currTrack));
                curvPeakNum(i,j) = numel(curvPeak);


                currYawNdx = yawPeakNdx(j);
                currTrack = ST{i}.all{j}(:,2);
                PeakYawPerc(i,j) = (currYawNdx-currTrack(1)+1)./length(currTrack);
                yawPeakDelay(i,j) = tmp(j);
                valYaw{i,j} = currYaw(currTrack);
                valCurv{i,j} = currCurv(currTrack).*180./pi;
                valDiffYawCurv{i,j} = valCurv{i,j}-valYaw{i,j};
                valSlipAng{i,j} = mvmtAng(currTrack);
                valDiffSlipAng{i,j} = dmvmtAng(currTrack);
                percBackwards(i,j) = sum(abs(mvmtAng(currTrack))>90)./length(currTrack);
                distance(i,j) = sqrt(diff(sAll{i}.Center.x(currTrack([1,end]))).^2+diff(sAll{i}.Center.y(currTrack([1,end]))).^2);
            end
        end
        features.PeakCurvPerc = PeakCurvPerc;
        features.lenTrack = lenTrack;
        features.spdTrack = spdTrack;
        features.curvPeakNum = curvPeakNum;

        features.PeakYawPerc = PeakYawPerc;
        features.yawPeakDelay = yawPeakDelay;
        features.valYaw = valYaw;
        features.valCurv = valCurv;
        features.valDiffYawCurv = valDiffYawCurv;
        features.valSlipAng = valSlipAng;
        features.valDiffSlipAng = valDiffSlipAng;
        features.percBackwards = percBackwards;
        features.distance = distance;
    end
end




