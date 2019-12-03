function [empFlys,synthFlys,successRate] = RunAndTumbleWrapper(inputFName,options,fileOpts)
% This function is the wrapper code to set up the parameters to run the run
% and tumble model
%
% Inputs:
%    inputFName: Data file with all of the empirical data.
%    options: various options (true = yes)
%       options.kin = consider kinematic changes?
%       options.stops = consider stops?
%       options.curvedWalks = consider walks to have curvature? false=all
%           walks have a curvature of 0
%       options.borderChoice = consider border choice?
%       options.SharpTurnBiase = consider turn bias for sharp turns?
%       options.CurvedWalkBiase = consider turn bias for curved walks?
%       options.noise = add gaussian noise N(0,1) to tracks?
%    fileOpts: optional options
%       opts.border = normalized radial position of the light border (0>=b<=1)
%       opts.dataConsFold = folder with the consolidated data files
%
% Outputs:
%    empFlys: structure with empirical data
%       empFlys.x: nxframes matrix of x-body positions
%       empFlys.y: nxframes matrix of y-body positions
%       empFlys.r: nxframes matrix of radial body positions
%       empFlys.thrust: nxframes matrix of thrust
%       empFlys.slip: nxframes matrix of slip
%       empFlys.xH: nxframes matrix of x-head positions
%       empFlys.yH: nxframes matrix of y-head positions
%       empFlys.rH: nxframes matrix of radial head positions
%       empFlys.firstEntry: nx1 array of first entry frame for each fly
%       empFlys.stopsBefore: 1xn cell array of "long" stops before
%       empFlys.stopsDuring: 1xn cell array of "long" stops during
%    synthFlys: structure with synthetic fly data (same fields as empFlys)
%    successRate: probability of flies that fit the goodfly criterium

close all
load(inputFName,'empFlys','fs','stopWalkPDFs','turnPDFs','pStop'...
    ,'turnInDuringProb','durProbDetrend','durProb','ratInOut','ratOutIn');
totFly = 500;
border = options.border;
excludeNonEntrys = true;

%--------------------------------------------------------------------------
% for nans in the turn in probability, set to the closest non-nan value
for i = 1:2
    tmp = find(isnan(turnInDuringProb(i,:)));
    [closestPtBef,closestPtAft,~] = findBeforeAfter(tmp,...
        find(~isnan(turnInDuringProb(i,:))),'both');
    turnInDuringProb(i,tmp) = turnInDuringProb(i,min([closestPtBef;closestPtAft]));
end
%--------------------------------------------------------------------------
% % set nans to 50/50
% turnInDuringProb(isnan(turnInDuringProb)) = 0.5;
% turnInDuringProb = turnInDuringProb(:,2:end);
% turnInDuringProb(turnInDuringProb==0) = 0.5;
%--------------------------------------------------------------------------

% setup the stopping probability
turnInProb = cell(1,8);turnInWalkProb = cell(1,8);
if options.stops == false
    pStop.before = 0;
    pStop.d_o = 0;
    pStop.d_i = 0;
end

% if not considering changes in kinematics, then use the before kinematics
% in the during inside and during outside scenarios
if options.kin == false
    for i = 2:3
        stopWalkPDFs.prob{i} = stopWalkPDFs.prob{1};
        stopWalkPDFs.val{i} = stopWalkPDFs.val{1};
    end
    for i = 5:6
        stopWalkPDFs.prob{i} = stopWalkPDFs.prob{4};
        stopWalkPDFs.val{i} = stopWalkPDFs.val{4};
    end
    stopWalkPDFs.prob{8} = stopWalkPDFs.prob{7};
    stopWalkPDFs.val{8} = stopWalkPDFs.val{7};
    
    for i = 2:3
        turnPDFs.prob{i} = turnPDFs.prob{1};
        turnPDFs.val{i} = turnPDFs.val{1};
    end
    for i = 5:6
        stopWalkPDFs.prob{i} = turnPDFs.prob{4};
        turnPDFs.val{i} = turnPDFs.val{4};
    end
end
% if no turn bias, then change the bias to 50/50, otherwise use the bias
if options.SharpTurnBiase == false
    for i = 1:8
        turnInProb{1,i} = 0.5.*ones(size(turnInDuringProb));
    end 
else
    for i = 1:8
        turnInProb{1,i} = turnInDuringProb;
    end
    turnInProb{1,1} = 0.5.*ones(size(turnInDuringProb));
    turnInProb{1,4} = 0.5.*ones(size(turnInDuringProb));
end
if options.CurvedWalkBiase == false
    for i = 1:8
        turnInWalkProb{1,i} = 0.5.*ones(size(turnInDuringProb));
    end 
else
    for i = 1:8
        turnInWalkProb{1,i} = turnInDuringProb;
    end
    turnInWalkProb{1,1} = 0.5.*ones(size(turnInDuringProb));
    turnInWalkProb{1,4} = 0.5.*ones(size(turnInDuringProb));
end

% Normalize the detrended probability of turning in, just in case it wasn't
% previously normalized. Note that index 13 is the boundary between inside
% and outside
a = durProbDetrend(:,1:end-1);
% to account for flies moving too quickly for a sharp excess probability.
baseline = 0.2;
a(:,1:13) = a(:,1:13)./sum(a(:,1:13),2);
a(:,14:end) = a(:,14:end)./sum(a(:,14:end),2);

% initialize basic parameters
params.posInit = [0,0];
params.len = 10799;
params.pStop = pStop;
params.totFly = totFly;
params.posInit = [0,0];
BorderChoice.softBound = a+baseline;
BorderChoice.PerTracks = [ratInOut,ratOutIn];
PDFs.stopWalkPDFs = stopWalkPDFs;
PDFs.turnPDFs = turnPDFs;
turnIn.turnInProb = turnInProb;
turnIn.turnInWalkProb = turnInWalkProb;

% generate synthetic flies
[synthFlys,type,nTurn] = RunAndTumbleFinal(PDFs,turnIn,BorderChoice,params,options);

% choose good flies (first entry within 80th percentile of empirical and
% distance of at least 1.1xr and 0.9xr during stimulus on)
if excludeNonEntrys == true
    Q1 = quantile(empFlys.firstEntry,0.15);Q3 = quantile(empFlys.firstEntry,0.80);
    [synthFlys,type,nTurn,successRate] = excludeFlies(synthFlys,type,nTurn,border,Q1,Q3);
end
% generate the four kinematic conditions for the synthetic flies
[curvPks,curvWalks,stopCond,boundCond] = genSTCWSynth2(synthFlys,type);

% plottingRadChoiceTurnBias(s,['Figures\' fName(1:find(fName=='\')-1) '\RadChoiceTurnBias'])
% load('BestFit.mat','GlobMinX');
% [curvPks,curvWalks,stopCond,boundCond,err] = genSTCWSynth(synthFlys,GlobMinX,fileOpts);

% save the data
save([options.fileName '.mat'],'synthFlys','empFlys','border','fs','type','nTurn',...
    'successRate','options','curvPks','curvWalks','stopCond','boundCond');

end

function [synthFlysNew,type,nTurn,successRate] = excludeFlies(synthFlys,type,nTurn,border,Q1,Q3)
%cutoff = 7200;
% good flies have to enter before Q3
goodFlies = synthFlys.firstEntry<Q3;% & synthFlys.firstEntry>5400;
goodFlies2 = false(size(synthFlys.firstEntry));
% good flies have to enter before Q3
for i = 1:length(goodFlies)
    goodFlies2(i,1) = sum(synthFlys.r(i,1:synthFlys.firstEntry(i))<border*0.9,2)>1 &...
    sum(synthFlys.r(i,synthFlys.firstEntry(i):end)<border*0.9,2)>1 &...
    sum(synthFlys.r(i,1:synthFlys.firstEntry(i))>border*1.1,2)>1 &...
    sum(synthFlys.r(i,synthFlys.firstEntry(i):end)>border*1.1,2)>1;
end
goodFlies = goodFlies & goodFlies2;

% update the synthetic flies structure to only keep goodflies
synthFlysNew.x = synthFlys.x(goodFlies,:);synthFlysNew.xH = synthFlys.x(goodFlies,:);
synthFlysNew.y = synthFlys.y(goodFlies,:);synthFlysNew.yH = synthFlys.y(goodFlies,:);
synthFlysNew.r = synthFlys.r(goodFlies,:);synthFlysNew.rH = synthFlys.r(goodFlies,:);
synthFlysNew.thrust = sqrt(diff(synthFlysNew.x,[],2).^2+diff(synthFlysNew.x,[],2).^2).*30.*4;
synthFlysNew.slip = zeros(size(synthFlysNew.r));
synthFlysNew.firstEntry = synthFlys.firstEntry(goodFlies,:);
synthFlysNew.stopsBefore = synthFlys.stopsBefore(goodFlies);
synthFlysNew.stopsDuring = synthFlys.stopsDuring(goodFlies);

type = type(goodFlies,:);
nTurn = nTurn(goodFlies,:);

% calculate the % of flies that are considered good
successRate = sum(goodFlies)./length(goodFlies);
end





