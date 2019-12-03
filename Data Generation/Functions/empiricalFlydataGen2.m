function [empFlys,stopWalkPDFs,turnPDFs,rawVals,pStop,gen] = empiricalFlydataGen2(file,opts)
% This function is a wrapper function that generates the probability
% distributions for each of the 4 states
%
% Inputs:
%    file1: consolidated data file with the following information
%       'synthFlys','empFlys','curvPks','curvWalks','stopCond','boundCond','gen'
%    opts: structure with field border for light border
% 
% Outputs:
%    empFlys: structure with fields: x,y,r,xH,yH,rH, and firstEntry 
%    stopWalkPDFs, turnPDFs have the following format:
%       pdfs: Structure of the pmfs for each scenario
%           pdfs.type: 1xn cell array to label what each dimension means
%           pdfs.prob: linearized probability density (m grids x n dimensions)
%           pdfs.val: corresponding values for each grid to pdfs.prob
%    rawVals: structure the empirical values that are used to generate the
%       pmfs
%    pStop: structure with probability of stopping before, during inside,
%       and during outside
%    gen: string stating the genotype of the flies

load(file,'synthFlys','empFlys','curvPks','curvWalks','stopCond','boundCond','gen')

if exist('synthFlys','var')
    empFlys = synthFlys;
end
if ~exist('gen','var')
    gen = [];
end

% extract the kinematic values of the 4 movement states for before, during 
% inside and during outside
[empFlys,spdDist,durDist,durDist2,yawDist,curvDist,pStop] = ...
    rawValExtraction(empFlys,curvPks,curvWalks,stopCond,boundCond,opts);

% generate the joint pmf for each of the four states for before, during
% inside, and during outside
stopWalkPDFs = DistDurCurvRelationship2(spdDist,durDist,curvDist);
turnPDFs = turnPDFGen(yawDist,durDist2);

rawVals.spd = spdDist;
rawVals.dur = durDist;
rawVals.sharpTurn = yawDist;
rawVals.curvWalk = curvDist;

end