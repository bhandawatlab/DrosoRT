function [xModel,yModel,ST,err2] = genAbstractionTrack(s,curv,Cstarts,Cstops)
% This function takes start and end positions of sharp turns and creates an
% model of the track using sharp turn and curved walk abstractions
%
% Inputs:
%    Cs: Cell array with s.x and s.y (x and y positions of each track)
%    curv: cell array with curvature of each corresponding track in s
%    Cstarts: vector of sharp turn start index
%    Cstops: vector of sharp turn end index
% 
% Outputs:
%    xModel: modeled x positions
%    yModel: modeled y positions
%    ST: Structure with information about each sharp turn
%    err2: error

len = length(curv);
spd = sqrt(diff(s.x).^2+diff(s.y).^2);
x = s.x;
y = s.y;

turnLen = Cstops-Cstarts;
Cstops(turnLen==0) = [];
Cstarts(turnLen==0) = [];

turns2 = Cstarts(2:end)-Cstops(1:end-1);
ndx = find(turns2<2);
Cstarts(ndx+1) = [];Cstops(ndx) = [];

% initialize start and end indexes of straight tracks
if ~isempty(Cstops)
    straightTracks = Cstops(2:end)-Cstarts(1:end-1);
    straightTracksEnd = Cstarts([false straightTracks>0]);
    straightTracksStart = Cstops([straightTracks>0 false]);
    
    straightTracksStart = straightTracksStart+1;
    straightTracksEnd = straightTracksEnd-1;
    
else
    straightTracksStart = 1;straightTracksEnd = len;
end

if min(Cstarts)~=1
    straightTracksStart = [1 straightTracksStart];
    straightTracksEnd = [Cstarts(1)-1 straightTracksEnd];
end
if max(Cstops)~=len
    straightTracksStart = [straightTracksStart Cstops(end)+1];
    straightTracksEnd = [straightTracksEnd len];
end

% compute the length and average change in curvature of sharp turns
turnLen = Cstops-Cstarts+1;
if ~isempty(turnLen)
    for i = 1:length(Cstarts)
        tmp = abs(curv(Cstarts(i):Cstops(i)));
        maxTurn(i) = max(tmp);
        CurvNdx(i) = Cstarts(i)+ceil(mean(find(tmp==maxTurn(i))))-1;
        totTurn(i) = sum(curv(Cstarts(i):Cstops(i)));
    end
    avgTurn = totTurn./turnLen;
else
    totTurn = [];avgTurn = [];CurvNdx = [];maxTurn = [];
end

% compute the track length and average change in curvature of curved
% (constant smooth curving) tracks
trackLen = straightTracksEnd-straightTracksStart+1;
for i = 1:length(straightTracksEnd)
    maxCurv(i) = sum(abs(curv(straightTracksStart(i):straightTracksEnd(i))));
    totCurv(i) = sum(curv(straightTracksStart(i):straightTracksEnd(i)));
end
avgCurv = totCurv./(trackLen);

% generate error and figures
ST.Curv = totTurn;
ST.Cstarts = Cstarts;
ST.Cstops = Cstops;
ST.len = turnLen;
ST.ndx = CurvNdx;
CW.Curv = avgCurv;
CW.Cstarts = straightTracksStart;
CW.Cstops = straightTracksEnd;
CW.len = trackLen;

[xModel,yModel,~] = generateSyntheticTrack(s,spd,curv,ST,CW);
err2 = (sqrt(mean((xModel-x).^2))+sqrt(mean((yModel-y).^2)));

end

function [xModel,yModel,synCurv] = generateSyntheticTrack(s,spd,curv,ST,CW)

% generate a sequence of synthetic curvature where the fly turns sharp
% turns in one frame and has a steady turn at straight tracks
lenST = length(ST.len);lenCW = length(CW.len);lenTrack = length(spd);
avgSpd = zeros(size(spd));synCurv = zeros(1,lenTrack);
for i = 1:lenST
    if ST.Cstops(i) == lenTrack
        tmpSeg = ST.Cstarts(i):ST.Cstops(i);
    else
        tmpSeg = ST.Cstarts(i):ST.Cstops(i);
    end
    avgSpd(tmpSeg) = avgSpd(tmpSeg)+mean(spd(tmpSeg));
end
synCurv(ST.ndx) = ST.Curv;

for i = 1:lenCW
    if CW.Cstops(i) == lenTrack
        tmpSeg = CW.Cstarts(i):CW.Cstops(i);
    else
        tmpSeg = CW.Cstarts(i):CW.Cstops(i);
    end
    avgSpd(tmpSeg) = mean(spd(tmpSeg));
    synCurv(tmpSeg) = synCurv(tmpSeg)+CW.Curv(i);
end


dXY = [diff(s.x);diff(s.y)];
Ang = myatan(dXY(1,:),dXY(2,:),'degrees',2);
Ang2 = Ang(1)+[0 (cumsum((curv(1:end-1)+curv(2:end)).*180./pi./2))];
Ang3 = Ang(1)+[0 (cumsum((synCurv(1:end-1)+synCurv(2:end)).*180./pi./2))];

xModel = s.x(1)+[0 cumsum(avgSpd.*cosd(Ang3))];
yModel = s.y(1)+[0 cumsum(avgSpd.*sind(Ang3))];

% xModel2 = s.x(1)+[0 cumsum(spd.*cosd(Ang))];
% yModel2 = s.y(1)+[0 cumsum(spd.*sind(Ang))];

end





















