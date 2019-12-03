function [empFlys,spdDist,durDist,durDist2,yawDist,curvDist,pStop] = ...
    rawValExtraction(empFlys,curvPks,curvWalks,stopCond,boundCond,opts)
% This function extracts the kinematic distributions of the 4 movement 
% states for before, during inside and during outside
%
% Inputs:
%    empFlys: structure with fields: x,y,r,xH,yH,rH, and firstEntry 
%    curvPks: structure detailing sharp turn condition, specifically how
%       long it occured, and curvature
%    curvWalks: structure detailing curved walk condition, specifically how
%       long it occured, and curvature
%    stopCond: structure detailing stop condition, specifically how long it
%       occured, and amount turned
%    boundCond: structure detailing boundary condition, specifically how
%       long it occured, and radial angle movement
%    opts: structure with field border for light border
% 
% Outputs:
%    empFlys: structure with fields: x,y,r,xH,yH,rH, and firstEntry 
%    spdDist: 1x8 cell matrix with the speed for each of the 8 scenarios
%    durDist: 1x8 cell matrix with the duration for each of the 8 scenarios
%       contains curved walks, but not sharp turns
%    durDist2: 1x8 cell matrix with the duration for each of the 8
%       scenarios contains sharp turns, but not curved walks
%    yawDist: 1x8 cell matrix with total curvature for sharp turns and
%       stops
%    curvDist: 1x8 cell matrix with total curvature for curved walks
%    pStop: structure with probability of stopping before, during inside,
%       and during outside

spdDist = cell(1,8);
durDist = cell(1,8);durDist2 = cell(1,8);
yawDist = cell(1,8);
curvDist = cell(1,8);

% define speed and duration of experiment
totFly = length(empFlys.firstEntry);
spdAll = sqrt(diff(empFlys.x.*10,[],2).^2+diff(empFlys.y.*10,[],2).^2)';
trialDur = size(spdAll,1);

% normalizes position based on radius of 4 cm
[empFlys] = normalizeDat(empFlys,4);

% find when the fly is inside and the light is turned on
in = empFlys.rH'<opts.border;
during = false(size(in));
for i = 1:totFly
    during(empFlys.firstEntry:end,i) = true; 
end

% define before, during out, and during in scenarios
allScenarios(:,:,1) = ~during;      % before
allScenarios(:,:,2) = during & ~in; % during out
allScenarios(:,:,3) = during & in;  % during in

% loop through each fly
for i = 1:totFly
    % sharp turns
    [yawDistTmp,durDistTmp] = getTurnDur(curvPks,allScenarios,i,trialDur);
    for j = 1:3
        yawDist{1,j+3} = [yawDist{1,j+3} mod(yawDistTmp{j},360)];
        durDist2{1,j+3} = [durDist2{1,j+3} durDistTmp{j}];
    end

    % stops
    [yawDistTmp,durDistTmp] = getTurnDur(stopCond,allScenarios,i,trialDur);
    for j = 1:3
        curvDist{1,j} = [curvDist{1,j} mod(yawDistTmp{j},360)];
        yawDist{1,j} = [yawDist{1,j} mod(yawDistTmp{j},360)];
        durDist{1,j} = [durDist{1,j} durDistTmp{j}];
        durDist2{1,j} = [durDist2{1,j} durDistTmp{j}];
    end
    
    % boundary
    [yawDistTmp,durDistTmp] = getTurnDur(boundCond,allScenarios,i,trialDur);
    for j = 1:2
        curvDist{1,j+6} = [curvDist{1,j+6} mod(yawDistTmp{j},360)];
        durDist{1,j+6} = [durDist{1,j+6} durDistTmp{j}];
    end
    
    % curved walks
    currAng = abs(curvWalks.avg{1,i});
    [currDur,~] = cellfun(@size,curvWalks.all{1,i});
    currSpd = zeros(1,length(currDur));Per = zeros(length(currDur),3);
    for j = 1:length(currDur)
        tmp = curvWalks.all{1,i}{1,j}(:,2);
        currSpd(j) = mean(spdAll(tmp(tmp<=trialDur),i));
        tmp2 = allScenarios(tmp(tmp<=trialDur),i,:);
        z = size(tmp2);
        Per(j,:) = sum(reshape(tmp2,[z(1) z(3)]))./currDur(j);
    end
    [~,ndx] = sort(Per,2,'descend');
    type = ndx(:,1);
    % assign curvWalk variables to the different scenarios
    for j = 1:3
        curvDist{1,j+3} = [curvDist{1,j+3} mod(currAng(type==j).*180./pi,360)];
        durDist{1,j+3} = [durDist{1,j+3} currDur(type==j)];
        spdDist{1,j+3} = [spdDist{1,j+3} currSpd(type==j)];
    end
end

for i = 1:3
    spdDist{1,i} = 0;
end

% probabilities of stopping
pStop.before = length(durDist{1})./(length(durDist{1})+length(durDist{4}));
pStop.d_o = length(durDist{2})./(length(durDist{2})+length(durDist{5}));
pStop.d_i = length(durDist{3})./(length(durDist{3})+length(durDist{6}));

end

function [yawDist,durDist] = getTurnDur(cond,allScenarios,flyN,trialDur)
yawDist = cell(1,3);durDist = cell(1,3);
if ~isempty(cond.all{1,flyN})
    % only consider scenarios that are greater tham 3 frames in length
    cond.all{1,flyN} = cond.all{1,flyN}(cellfun(@length,cond.all{1,flyN})>3);
    cond.ndx{1,flyN} = cond.ndx{1,flyN}(cellfun(@length,cond.all{1,flyN})>3);
    
    [currDur,~] = cellfun(@size,cond.all{1,flyN});
    currAng = abs(cond.tot{1,flyN}(cond.ndx{1,flyN}<=trialDur));
    tmp = cond.ndx{1,flyN}(cond.ndx{1,flyN}<=trialDur);
    tmp2 = allScenarios(tmp,flyN,:);
    z = size(tmp2);
    [~,type] = sort(reshape(tmp2,[z(1) z(3)]),2,'descend');
    type = type(:,1);
    for j = 1:3
        yawDist{1,j} = currAng(type==j).*180./pi;
        durDist{1,j} = currDur(type==j);
    end
end

end





