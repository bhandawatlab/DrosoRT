function [SharpCurv,CurvTrack,err2] = modelCurvedWalks(s,curv,Cstarts,Cstops,pltInf)
% This function takes start and end positions of sharp turns and creates an
% model of the track using sharp turn and curved walk abstractions
%
% Inputs:
%    s: Cell array with s.x and s.y (x and y positions of each track)
%    curv: cell array with curvature of each corresponding track in s
%    Cstarts: vector of sharp turn start index
%    Cstops: vector of sharp turn end index
%    pltInf: subplot information [x,y,n]. No plotting if pltInf=[]
% 
% Outputs:
%    SharpCurv: Structure with information about each sharp turn
%    CurvTrack: Structure with information about each curved walk
%    err2: error

% get some basic data
len = length(curv);
spd = sqrt(diff(s.x).^2+diff(s.y).^2);
x = s.x;
y = s.y;

% any turns that has a duration of 1 frame is deleted
turnLen = Cstops-Cstarts;
Cstops(turnLen==0) = [];
Cstarts(turnLen==0) = [];
% any turns that has a duration less than 3 frame is deleted (this is up
% for debate)
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

% define some structures to organize data
ST.Curv = totTurn;
ST.Cstarts = Cstarts;
ST.Cstops = Cstops;
ST.len = turnLen;
ST.ndx = CurvNdx;
CW.Curv = avgCurv;
CW.Cstarts = straightTracksStart;
CW.Cstops = straightTracksEnd;
CW.len = trackLen;

% sum(totTurn)+sum(avgCurv.*trackLen);
% sum(curv)

% generate modeled tracks 
[xModel,yModel,~] = generateSyntheticTrack(s,spd,curv,ST,CW,pltInf);
% calculate error
err2 = (sqrt(mean((xModel-x).^2))+sqrt(mean((yModel-y).^2)));

% define sharp turn matrices
SharpCurv.tot = totTurn;
SharpCurv.avg = avgTurn;
SharpCurv.dur = turnLen;
SharpCurv.max = maxTurn;
SharpCurv.ndx = CurvNdx;
SharpCurv.len = turnLen;
if ~isempty(turnLen)
    for i = 1:length(Cstarts)
        SharpCurv.all{i} = [curv(Cstarts(i):Cstops(i));Cstarts(i):Cstops(i)]';
    end
    % compute angle entering trajectory and angle leaving trajectory
    [SharpCurv] = computeTurnDirections(x,y,Cstarts,Cstops,CurvNdx,SharpCurv);
else
    SharpCurv.all = {};SharpCurv.dirRelCenterBef = [];SharpCurv.dirRelCenterAft = [];
end


% define curved walk matrices
CurvTrack.tot = totCurv;
CurvTrack.avg = avgCurv;
CurvTrack.dur = trackLen;
CurvTrack.max = maxCurv;
CurvTrack.ndx = round((straightTracksEnd+straightTracksStart)./2);
CurvTrack.len = trackLen;
for i = 1:length(straightTracksEnd)
    CurvTrack.all{i} = [curv(straightTracksStart(i):straightTracksEnd(i));
        straightTracksStart(i):straightTracksEnd(i)]';
end
% compute angle entering trajectory and angle leaving trajectory
[CurvTrack] = computeTurnDirections(x,y,straightTracksStart,straightTracksEnd,CurvTrack.ndx,CurvTrack);


end

function [xModel,yModel,synCurv] = generateSyntheticTrack(s,spd,curv,ST,CW,pltInf)

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

% calculate the xy and movement angle of the tracks
dXY = [diff(s.x);diff(s.y)];
Ang = myatan(dXY(1,:),dXY(2,:),'degrees',2);
Ang2 = Ang(1)+[0 (cumsum((curv(1:end-1)+curv(2:end)).*180./pi./2))];
Ang3 = Ang(1)+[0 (cumsum((synCurv(1:end-1)+synCurv(2:end)).*180./pi./2))];

% x = s.x(1)+[0 cumsum(spd.*cosd(Ang2))];
% y = s.y(1)+[0 cumsum(spd.*sind(Ang2))];

% calculate the modeled tracks
xModel = s.x(1)+[0 cumsum(avgSpd.*cosd(Ang3))];
yModel = s.y(1)+[0 cumsum(avgSpd.*sind(Ang3))];

% plot the empirical and modeled tracks
if ~isempty(pltInf)
    xmin = min([s.x xModel -4]);ymin = min([s.y yModel -4]);
    xmax = max([s.x xModel 4]);ymax = max([s.y yModel 4]);
    
    subplot(pltInf(1),pltInf(2),pltInf(3));
    plot(s.x,s.y,'k','LineWidth',2);hold on
    for i = 1:length(ST.Cstarts)
        plot(s.x(ST.Cstarts(i):ST.Cstops(i)),s.y(ST.Cstarts(i):ST.Cstops(i)),'r','LineWidth',2)
    end
    axis([xmin xmax ymin ymax])
    subplot(pltInf(1),pltInf(2),pltInf(3)+1);hold on
    plot(xModel,yModel,'k','LineWidth',2)
    axis([xmin xmax ymin ymax])
end

end





















