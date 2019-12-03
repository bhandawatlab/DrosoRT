function [errAll,SharpCurv,CurvTrack] = objFunc(p1,p2,p3,p4,p5,p6,s,curv,pltInf)
% This is the objective function to calculate the parameters needed to
% delineate between sharp turns and curved walks
%
% Inputs:
%    p1: Curvature Weight
%    p2: Change in Curvature Weight
%    p3: Curvature Threshold (in radians)
%    p4: Change in Curvature Threshold (in radians)
%    p5: Sharp turn threshold multiple
%    p6: Curved walk threshold multiple
%    s: Cell array with s.x and s.y (x and y positions of each track)
%    curv: cell array with curvature of each corresponding track in s
%    pltInf: subplot information [x,y,n]. No plotting if pltInf=[]
% 
% Outputs:
%    errAll: total error
%    SharpCurv: Structure with information about each sharp turn
%    CurvTrack: Structure with information about each curved walk

param = [p1,p2,p3,p4,p5,p6];
%tic
l = length(curv);distance = zeros(1,l);
Cstarts = cell(1,l);Cstops = cell(1,l);errPos = zeros(1,l);
% loop through each track
for c = 1:l
    % calculate where each sharp turn/ curved walk trajectory starts/ends 
    % based on the current parameters
    [Cstarts{c},Cstops{c},~] = FindCurvStartStopP(curv{c},param);
    % short trajectories are removed
    shortST = Cstops{c}-Cstarts{c}<2;
    Cstarts{c}(shortST) = [];Cstops{c}(shortST) = [];
    if length(curv{c})>3
        % create modeled tracks based on the sharp turn and curved walk
        % start/end positions calculated above and then calculate the error
        [SharpCurv,CurvTrack,errPos(c)] = modelCurvedWalks(s{c},curv{c},Cstarts{c},Cstops{c},pltInf);
    else
        [SharpCurv,CurvTrack] = initiateEmpty();
    end
    % calculate distance
    distance(c) = sum(sqrt(diff(s{c}.x).^2+diff(s{c}.y).^2));
end
% get average percentage error normalized to the distance
errAll = sum(errPos./distance)./l.*100;% avg percent error
%errAll = sum(errPos);
%toc
end

function [SharpCurv,CurvTrack] = initiateEmpty()

% define sharp turn matrices
SharpCurv.tot = [];
SharpCurv.avg = [];
SharpCurv.dur = [];
SharpCurv.max = [];
SharpCurv.ndx = [];
SharpCurv.all = {};
SharpCurv.dirRelCenterBef = [];
SharpCurv.dirRelCenterAft = [];

% define curved walk matrices
CurvTrack.tot = [];
CurvTrack.avg = [];
CurvTrack.dur = [];
CurvTrack.max = [];
CurvTrack.ndx = [];
CurvTrack.all = {};
CurvTrack.dirRelCenterBef = [];
CurvTrack.dirRelCenterAft = [];
end