function [features,curvPks,curvWalks] = processFeaturesKatsov(yawThresh,curvThresh,fitFile)
% This function generates features of each sharp turn such as yaw, curv,
% etc from the Katsov data set L and thresholds for yaw and curvature
% 
% Inputs:
%    yawThresh: threshold for finding peaks in yaw
%    curvThresh: threshold for finding peaks in curvatures
%    fitFile: file with the best fit parameters to delineate between sharp
%       turns and curved walks
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
%    curvPks: struction with information about sharp turns
%    curvWalks: struction with information about curved walks

load('Data Full/KatsovData/dataset_L.mat','UI','VR','VRTH','VS','VT')
fs = 30;
load(fitFile,'GlobMinX');
err = zeros(1,max(UI));
figNum = 1;

nCt = max(UI);maxL = 0;
curvPks = cell(1,nCt);curvWalks = cell(1,nCt);
YawPks = cell(1,nCt);curvAll = cell(1,nCt);
mvmtAngAll = cell(1,nCt);dmvmtAngAll = cell(1,nCt);
xAll = cell(1,nCt);yAll = cell(1,nCt);
progressbar
for i = 1:nCt
    currTrack = find(UI==i);
    
    thrust = VT(currTrack)';
    slip = VS(currTrack)';
    yaw = VR(currTrack)'./fs;
    
    initAng = myatan(thrust(1),slip(1),'degrees',2);
    initAng = initAng+VRTH(currTrack(1))./30;
    initX = 0;
    initY = 0;
    
    allAng = initAng+[0 cumsum(yaw)];
    allAng = allAng-360*floor(allAng./360);
    
    xThrust = cosd(allAng(1:end-1)).*thrust./fs;
    yThrust = sind(allAng(1:end-1)).*thrust./fs;
    
    xSlip = cosd(allAng(1:end-1)-90).*slip./fs;
    ySlip = sind(allAng(1:end-1)-90).*slip./fs;
    
    delX = xThrust+xSlip;
    delY = yThrust+ySlip;
    
    xRec = initX+[0 cumsum(delX)];
    yRec = initY+[0 cumsum(delY)];
    
    mvmtAng = allAng(2:end)-myatan(diff(xRec),diff(yRec),'degrees',2);
    mvmtAng(mvmtAng>180) = mvmtAng(mvmtAng>180)-floor(mvmtAng(mvmtAng>180)./180).*180;
    mvmtAng(mvmtAng<-180) = mvmtAng(mvmtAng<-180)-ceil(mvmtAng(mvmtAng<-180)./180).*180;
    dmvmtAng = [0 diff(mvmtAng)];
    
    mvmtAngAll{i} = mvmtAng;
    dmvmtAngAll{i} = dmvmtAng;
    xAll{i} = xRec;yAll{i} = xRec;
    
    % find peaks in yaw
    [~,pks] = findpeaks(abs(yaw).*pi./180,'MinPeakHeight',yawThresh);
    YawPks{i} = pks;
    
    [curv,~] = CalcCurvatureKatsov(xRec,yRec);
    curvAll{i} = curv;
    
    s.x = xRec;
    s.y = yRec;
    
    if mod(i,10) == 1
        %figNum = figNum+1;figure(figNum);
        %set(gcf,'Position',[2 42 638 600])
        k = 1;
    end
    
    pltInf = [];%[5,4,k];
    obj = @(x)objFunc(x(1),x(2),x(3),x(4),x(5),x(6),{s},{curv},pltInf);
    [err(i),curvPksTmp,curvWalksTmp] = obj(GlobMinX);
    
    if ~isempty(curvPksTmp) && ~isempty(curvWalksTmp)
        curvPks{i} = curvPksTmp;
        curvWalks{i} = curvWalksTmp;
        maxL = max(maxL,length(curvPks{i}.ndx));
    end
    k = k+2;
    progressbar(i/nCt)
end


PeakCurvPerc = zeros(nCt,maxL);lenTrack = zeros(nCt,maxL);
spdTrack = zeros(nCt,maxL);curvPeakNum = zeros(nCt,maxL);
PeakYawPerc = zeros(nCt,maxL);yawPeakDelay = zeros(nCt,maxL);
percBackwards = zeros(nCt,maxL);distance = zeros(nCt,maxL);
valYaw = cell(nCt,maxL);valSlip = cell(nCt,maxL);
valCurv = cell(nCt,maxL);valSlipAng = cell(nCt,maxL);
valDiffSlipAng = cell(nCt,maxL);valDiffYawCurv = cell(nCt,maxL);
progressbar
for i = 1:nCt
    if ~isempty(curvPks{i})
        currTrack = find(UI==i);

        % get kinematics
        currYaw = VR(currTrack)'./fs;
        currSlip = VS(currTrack)'./fs;
        currThrust = VT(currTrack)'./fs;
        currCurv = curvAll{i};
        spd = sqrt(currThrust.^2+currSlip.^2);
        
        mvmtAng = mvmtAngAll{i};
        dmvmtAng = dmvmtAngAll{i};

        % obtain delays between yaw peak and the sharp turn
        currYawPks = YawPks{i};
        currNdx = curvPks{i}.ndx;
        
        if isempty(currYawPks)
            tmp = 1000.*ones(size(currNdx));
            yawPeakNdx = currNdx+tmp;
        else
            [closestPtBef,closestPtAft,~] = findBeforeAfter(currNdx,currYawPks,'both');
            befDiff = closestPtBef-currNdx;
            aftDiff = closestPtAft-currNdx;
            [~,I] = min([abs(aftDiff);abs(befDiff)],[],1);
            tmp = aftDiff;
            tmp(I==2) = befDiff(I==2);
            [~,ia,~] = intersect(currNdx,currYawPks);
            tmp(ia) = 0;
            yawPeakNdx = currNdx+tmp;
        end


        for j = 1:length(curvPks{i}.ndx)
            currNdx = curvPks{i}.ndx(j);
            currTrack = curvPks{i}.all{j}(:,2);

            if length(currTrack)>2
                curvTmp = abs(curvAll{i}(currTrack));
                [~,curvPeak] = findpeaks(curvTmp,'MinPeakHeight',curvThresh);

                PeakCurvPerc(i,j) = (currNdx-currTrack(1)+1)./length(currTrack);
                lenTrack(i,j) = length(currTrack);
                spdTrack(i,j) = mean(spd(currTrack));
                curvPeakNum(i,j) = numel(curvPeak);
                
                currYawNdx = yawPeakNdx(j);
                currTrack = curvPks{i}.all{j}(:,2);
                PeakYawPerc(i,j) = (currYawNdx-currTrack(1)+1)./length(currTrack);
                yawPeakDelay(i,j) = tmp(j);
                valYaw{i,j} = currYaw(currTrack);
                valSlip{i,j} = currSlip(currTrack);
                valCurv{i,j} = currCurv(currTrack).*180./pi;
                valDiffYawCurv{i,j} = valCurv{i,j}-valYaw{i,j};
                valSlipAng{i,j} = mvmtAng(currTrack);
                valDiffSlipAng{i,j} = dmvmtAng(currTrack);
                percBackwards(i,j) = sum(abs(mvmtAng(currTrack))>90)./length(currTrack);
                distance(i,j) = sqrt(diff(xAll{i}(currTrack([1,end]))).^2+diff(yAll{i}(currTrack([1,end]))).^2);
            end
        end
    end
    progressbar(i/nCt)
end

features.PeakCurvPerc = PeakCurvPerc;
features.lenTrack = lenTrack;
features.spdTrack = spdTrack;
features.curvPeakNum = curvPeakNum;

features.PeakYawPerc = PeakYawPerc;
features.yawPeakDelay = yawPeakDelay;
features.valYaw = valYaw;
features.valSlip = valSlip;
features.valCurv = valCurv;
features.valDiffYawCurv = valDiffYawCurv;
features.valSlipAng = valSlipAng;
features.valDiffSlipAng = valDiffSlipAng;
features.percBackwards = percBackwards;
features.distance = distance;

end




