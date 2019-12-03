function [out] = computeStopBorder(vec,s,curv,ndx,out)
% This function creates structures for the stop or border state
%
% Inputs:
%    vec: vector of ones/true and zeros/false indicating whether or not in
%       that state
%    s: structure with fields Center.x, Center.y
%    curv: vector same length as vec that has the curvature values
%    ndx: current fly index
%    out: structure for the stop or border state with fields max, tot, dur,
%       avg, ndx, and all
% 
% Outputs:
%    out: structure for the stop or border state with fields max, tot, dur,
%       avg, ndx, and all

% find where the vector starts and ends for continuous segments of numbers
[startNdx,endNdx,type] = startEndSeq(vec);
startNdx = startNdx(type == 1);
endNdx = endNdx(type == 1);

x = s.Center.x;
y = s.Center.y;

% if no segment of ones, then initialize with empty array
if isempty(startNdx)
    out.max{ndx} = [];
    out.tot{ndx} = [];
    out.dur{ndx} = [];
    out.avg{ndx} = [];
    out.ndx{ndx} = [];
    out.all{ndx} = [];
    out.dirRelCenterBef{ndx} = [];
    out.dirRelCenterAft{ndx} = [];
else
    % otherwise loop through each segment
    for i = 1:length(startNdx)
        curvTmp = curv(startNdx(i):endNdx(i));
        out.max{ndx}(i) = max(curvTmp);
        out.tot{ndx}(i) = sum(curvTmp);
        out.dur{ndx}(i) = length(curvTmp);
        out.avg{ndx}(i) = mean(curvTmp);
        out.ndx{ndx}(i) = floor((startNdx(i)+endNdx(i))./2);
        out.all{ndx}{i} = [curvTmp',[startNdx(i):endNdx(i)]'];
        
        u = [x(out.ndx{ndx}(i))-x(startNdx(i));y(out.ndx{ndx}(i))-y(startNdx(i));zeros(1)];
        v = [x(endNdx(i))-x(out.ndx{ndx}(i));y(endNdx(i))-y(out.ndx{ndx}(i));zeros(1)];
        w = [-x(out.ndx{ndx}(i));-y(out.ndx{ndx}(i));zeros(1)];
        out.dirRelCenterBef{ndx}(i) = atan2(norm(cross(u,w)),dot(u,w)).*180./pi;
        out.dirRelCenterAft{ndx}(i) = atan2(norm(cross(v,w)),dot(v,w)).*180./pi;
    end
end
end