function [velClass2,velBinary2,runstops,genTypeNdx] = Velocity_Classifier...
    (velocity,stopthreshold,runthreshold)
% This function classifies the movements into moving and not moving based
% on stop and run thresholds
%
% Inputs:
%    velocity: speed of flies 1xg genotype, each cell is (m time points x n
%       flies)
%    stopthreshold: threshold for stopped
%    runthreshold: threshold for running
% 
% Outputs:
%    velClass2: same as the velocity input, but al instances where flys are
%       deemed stopped, the velocity is set to 0
%    velBinary2: an array of 1's for moving and 0's for stopped
%    runstops: transition points from moving to stops
%    genTypeNdx: indexing by genotypes

% define number of genotypes
numGenotype = length(velocity);

genTypeNdx = cell(1,length(velocity));
ndx = 1;
%define indexes (rows) for genotypes after consolidating velocity data to
%one matrix
for i = 1:numGenotype
    genTypeNdx{i} = ndx:ndx+size(velocity{i},1)-1;
    ndx = ndx+size(velocity{i},1);   
end
velAll = cell2mat(velocity');

%initialize matrices
velClass = velAll;
velBinary = nan(size(velAll));
runstops = zeros(size(velAll));

% declare first points for each fly as stop or run based on middle of
% stop/run threshold
velBinary(velAll(:,1) > (stopthreshold+runthreshold)./2,1) = 1;
velBinary(velAll(:,1) <= (stopthreshold+runthreshold)./2,1) = 0;


for i = 2:size(velAll,2)-1
    prevBin = velBinary(:,i-1);
    currVel = velAll(:,i);
    nextVel = velAll(:,i+1);
    
    case1 = (prevBin == 1) & (currVel >= stopthreshold);
    case2 = (prevBin == 1) & (currVel < stopthreshold) & (nextVel <= stopthreshold);
    case3 = (prevBin == 1) & (currVel < stopthreshold) & (nextVel > stopthreshold);
    
    case4 = (prevBin == 0) & (currVel < runthreshold);
    case5 = (prevBin == 0) & (currVel >= runthreshold);
    
    velBinary(case1,i) = 1; % previous frame = run, and current frame above stop threshold
    velBinary(case2,i) = 0; % previous frame = run, and current frame below stop threshold and next frame below stop threshold
    velBinary(case3,i) = 1; % previous frame = run, and current frame below stop threshold and next frame above stop threshold
    velBinary(case4,i) = 0; % previous frame = stop, and current frame below run threshold
    velBinary(case5,i) = 1; % previous frame = stop, and current frame above run threshold
    
    runstops(case2,i) = 1;  % transition from run to stop
    runstops(case5,i) = 1;  % transition from stop to run
end
% last frame for each track is the same as the second to last frame
velBinary(:,end) = velBinary(:,end-1);
velClass(velBinary==0) = 0;

% restructure data into cells by genotype
velBinary2 = cell(1,numGenotype);velClass2 = cell(1,numGenotype);
for i = 1:numGenotype
    velBinary2{i} = velBinary(genTypeNdx{i},:);
    velClass2{i} = velClass(genTypeNdx{i},:);
end


end


