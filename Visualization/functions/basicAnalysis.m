function [s] = basicAnalysis(Flys,totFly,border)
% This function computes some basic analysis such as radial probability,
% attraction index, time to return, time inside before leaving
%
% Inputs:
%    Flys: structure with fields
%       Flys.r: radial location
%       Flys.firstEntry: time of first entry
%    totFly: total number of flies
%    border: radial location of light border
%
% Outputs:
%    s: structure with fields
%       s.Prob: radial occupancy
%       s.timeToReturn.before: time to return before first entry
%       s.timeInToTransit.before: time inside before leaving before first entry
%       s.timeToReturn.during: time to return before first entry
%       s.timeInToTransit.during: time inside before leaving after first entry
%       s.attractNdx: proportion of time spent inside

rMax = max(Flys.r(:));
Flys.r = Flys.r./rMax;

firstEntry = Flys.firstEntry;
% radial probability
Before = [];During = [];
for i = 1:totFly
    Before = [Before Flys.r(i,1:firstEntry(i))];
    During = [During Flys.r(i,firstEntry(i):end)];
end
edges = 0:0.05:1;
[N,~] = histcounts(Before,edges);
N = N+eps;
Prob.before.y = N./sum(N);
Prob.before.x = edges(2:end)-0.025;
[N,~] = histcounts(During,edges);
N = N+eps;
Prob.during.y = N./sum(N);
Prob.during.x = edges(2:end)-0.025;

% attraction index
for i = 1:totFly
    totInside = sum(Flys.r(i,firstEntry(i):end)<border);
    totOutside = sum(Flys.r(i,firstEntry(i):end)>=border);
    attractNdx.during(i) = totInside./(totInside+totOutside);
    
    totInside = sum(Flys.r(i,1:firstEntry(i))<border);
    totOutside = sum(Flys.r(i,1:firstEntry(i))>=border);
    attractNdx.before(i) = totInside./(totInside+totOutside);
end

% time to return and time before leaving
timeToReturnAfter = cell(totFly,1);
timeInToTransitAfter = cell(totFly,1);
for i = 1:totFly
    inside = Flys.r(i,firstEntry(i):end)<border;
    leaveNdx = find(diff(inside)<0);
    if inside(1) == true
        enterNdx = [1 find(diff(inside)>0)];
    else
        enterNdx = [find(diff(inside)>0)];
    end
    if length(leaveNdx)>1
        timeToReturnAfter{i} = enterNdx(2:end)-leaveNdx(1:length(enterNdx)-1);
    end
    timeInToTransitAfter{i} = leaveNdx-enterNdx(1:length(leaveNdx));
end

timeToReturnBefore = cell(totFly,1);
timeInToTransitBefore = cell(totFly,1);
for i = 1:totFly
    inside = [1 Flys.r(i,1:firstEntry(i))<border];
    leaveNdx = find(diff(inside)<0);
    if inside(1) == true
        enterNdx = [1 find(diff(inside)>0)];
    else
        enterNdx = [find(diff(inside)>0)];
    end
    if length(leaveNdx)>1
        timeToReturnBefore{i} = enterNdx(2:end)-leaveNdx(1:min(length(enterNdx),length(enterNdx))-1);
    end
    timeInToTransitBefore{i} = leaveNdx-enterNdx(1:length(leaveNdx));
    if ~isempty(leaveNdx) && ~isempty(enterNdx) && leaveNdx(1) == enterNdx(1)
        timeInToTransitBefore{i} = timeInToTransitBefore{i}(2:end);
    end
end

s.Prob = Prob;
s.timeToReturn.before = timeToReturnBefore;
s.timeInToTransit.before = timeInToTransitBefore;
s.timeToReturn.during = timeToReturnAfter;
s.timeInToTransit.during = timeInToTransitAfter;
s.attractNdx = attractNdx;

end