function [synthFlys,state,nTurn] = RunAndTumbleFinal(PDFs,turnInProb,BorderChoice,params,options)
% This function is the code to run the run and tumble models
%
% Inputs:
%    PDFs: Data file with all of the empirical data.
%       stopWalkPDFs = structure for pdf of stop and curved walks tracks
%       turnPDFs = structure for pdf of sharp turns
%    turnInProb: various options (true = yes)
%       turnInProb = sharp turn turn in probability
%       turnInWalkProb = curved walk turn in probability
%    BorderChoice: optional options
%       softBound = probability of turning in
%       PerTracks = percent of tracks that had at least n turns before
%       hitting wall (nx2 matrix)
%    params: optional options
%       posInit = initial position that synthetic flies start from
%       len = number of frames for each synthetic flies
%       pStop = probability of stopping
%       totFly = number of synthetic flies to make
%    options: various options regarding the number of complexities of the
%       model (see RunAndTumbleWrapper.m)
%
% Outputs:
%    synthFlys: structure with empirical data
%       synthFlys.x: nxframes matrix of x-body positions
%       synthFlys.y: nxframes matrix of y-body positions
%       synthFlys.r: nxframes matrix of radial body positions
%       synthFlys.thrust: nxframes matrix of thrust
%       synthFlys.slip: nxframes matrix of slip
%       synthFlys.xH: nxframes matrix of x-head positions
%       synthFlys.yH: nxframes matrix of y-head positions
%       synthFlys.rH: nxframes matrix of radial head positions
%       synthFlys.firstEntry: nx1 array of first entry frame for each fly
%       synthFlys.stopsBefore: 1xn cell array of "long" stops before
%       synthFlys.stopsDuring: 1xn cell array of "long" stops during
%    state: what state is the fly in at each time point (stop, sharp turn, etc)
%    nTurn: number of turns since last exit for each fly at each time point


curvedRun = options.curvedWalks;
noise = options.noise;
borderChoice = options.borderChoice;
len = params.len;
border = options.border;
pStop = params.pStop;
totFly = params.totFly;
posInit = params.posInit;
softBound = BorderChoice.softBound;
PerTracks = [BorderChoice.PerTracks;zeros(300,2)];

% convert PDFs format
PDFs.all.type = [PDFs.turnPDFs.type(1:6) PDFs.stopWalkPDFs.type(4:8)];
PDFs.all.prob = [PDFs.turnPDFs.prob(1:6) PDFs.stopWalkPDFs.prob(4:8)];
PDFs.all.val = [PDFs.turnPDFs.val(1:6) PDFs.stopWalkPDFs.val(4:8)];
turnInProb.all = [turnInProb.turnInProb(1:6) turnInProb.turnInWalkProb(4:6)];

%assume that the fly always start off either in or near the center
xPos = zeros(totFly,len);
yPos = zeros(totFly,len);
rPos = zeros(totFly,len);
if ~isempty(posInit)
    xPos(:,1) = params.posInit(:,1);
    yPos(:,1) = params.posInit(:,2);
end
rPos(:,1) = sqrt(xPos(:,1).^2+yPos(:,1).^2);

% initialize parameters
currAng = zeros(totFly,1); allAng = zeros(totFly,len);
during = false(totFly,1);
currDur = zeros(totFly,1);
oldDur = zeros(totFly,1);
currCurv = zeros(totFly,1);
loc = zeros(totFly,len);
state = zeros(totFly,len);
currSpd = zeros(totFly,1);
currPhi = zeros(totFly,1);
newTrack = true(totFly,1);
duringAll = false(totFly,len);
stopLength = 90;
turn = false(totFly,1);
sBLim = (0:0.1:size(softBound,2)*0.1)./4;
sBLim(end) = sBLim(end)-0.0051;

newTrackState = 3*ones(totFly,1);% start in the curve walk state
state(:,1) = newTrackState;
nTurn = zeros(totFly,len);

turningBias = zeros(totFly,1);
dDur = zeros(1,9);
for i = 1:9
    tmp = diff(unique(PDFs.all.val{i}(:,1)));
    dDur(i) = tmp(1);
end

% generate synthetic tracks based on run and tumble model
progressbar
for i = 1:len-1
    bound = rPos(:,i)>0.985;
    inside = rPos(:,i)<border;
    outside = ~bound & ~inside;
    loc(bound,i) = 1;loc(inside,i) = 2;loc(outside,i) = 3;
    
    % this section is for the turn counter
    if i>1
        nTurn(:,i) = nTurn(:,i-1);              % initialize current turn counter as the previous turn number
        tmp = (loc(:,i)==2 & loc(:,i-1)==3) | (loc(:,i)==3 & loc(:,i-1)==2);
        nTurn(tmp & during,i) = 0;              % reset the turn counter if the fly crosses the light border in the during case
    end
    
    %----------------------------------------------------------------------
    % this section makes the fly choose a new track if it reaches the edge
    if i>1
        newLoc = (loc(:,i-1) ~= loc(:,i));
        newTrack(newLoc & loc(:,i)==1) = true;
    end
    
    % set flies to the during case after first entry
    if ~all(during) && i>len/2
        during(during==false) = inside(during==false);
        duringAll(:,i) = during;
    end
    before = ~during;
    
    % this section makes the fly choose a new track if the fly is in the
    % odor ring when the light suddenly turns on
    if i>1
        newLoc = (duringAll(:,i-1) ~= duringAll(:,i));
        newTrack(newLoc & loc(:,i)==2) = true;
    end
    
    %----------------------------------------------------------------------
    % this section makes the fly have an increase chance of choosing a new 
    % track right when entering/exiting the odor zone
    if borderChoice && sum(during)>0
        if sum((rPos(:,i)<sBLim(end) & rPos(:,i)>=sBLim(1)))>0
            pSamp = rand(totFly,1);
            
            v1 = [cosd(currAng),sind(currAng),zeros(totFly,1)];
            w = [-xPos(:,i),-yPos(:,i),zeros(totFly,1)];
            a = cross(v1,w);
            b = dot(v1',w')';
            c = sqrt(sum(a.^2,2));
            facingOut = atan2d(c,b)>90;
            facingOut = (facingOut & outside) | inside;
            
            for j = 1:size(softBound,2)
                bordLoc = (rPos(:,i)<sBLim(j+1) & rPos(:,i)>=sBLim(j));
                if sBLim(j)>=border
                    SampChoice = pSamp<=softBound(1,j).*PerTracks(nTurn(:,i)+1,1);
                    SampChoice2 = pSamp<=softBound(2,j).*PerTracks(nTurn(:,i)+1,1);
                else
                    SampChoice = pSamp<=softBound(1,j).*PerTracks(nTurn(:,i)+1,2);
                    SampChoice2 = pSamp<=softBound(2,j).*PerTracks(nTurn(:,i)+1,2);
                end
                newTrack(bordLoc & during & SampChoice & state(:,i-1)==3 & nTurn(:,i)<3) = true;
                newTrack(bordLoc & during & SampChoice2 & state(:,i-1)==3 & nTurn(:,i)>=3) = true;
                %newTrack(bordLoc & during & SampChoice & state(:,i-1)==3 & nTurn(:,i)<3 & facingOut) = true;% facingOut
                %newTrack(bordLoc & during & SampChoice2 & state(:,i-1)==3 & nTurn(:,i)>=3) = true;
            end
        end
        turningBias(during & nTurn(:,i)<3) = 1;
        turningBias(during & nTurn(:,i)>=3) = 2;
    end
    
    % if the current state ends, then a new track begins
    newTrack(currDur<=0) = true;
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % this section is for choosing what new state to enter
    % stop: 1, turn: 2, run: 3, boundary: 4
    %
    if i>1
        temp = rand(totFly,1);
        b = zeros(totFly,1);d_Inside = zeros(totFly,1);d_Outside = zeros(totFly,1);
        for s = 1:4
            b(state(:,i-1)==s & newTrack & before,1) = s;                         % before
            d_Inside(state(:,i-1)==s & newTrack & during & inside,1) = s;         % during inside
            d_Outside(state(:,i-1)==s & newTrack & during & outside,1) = s;       % during outside
        end

        % stops and sharp turns always go to runs
        for s = 1:2
            state(b==s, i) = 3;
            state(d_Inside==s, i) = 3;
            state(d_Outside==s, i) = 3;
        end

        % runs can go to stops
        state(temp<pStop.before & b==3, i) = 1;
        state(temp<pStop.d_o & d_Inside==3, i) = 1;
        state(temp<pStop.d_i & d_Outside==3, i) = 1;

        % runs that do not go to stops go to sharp turns
        state(temp>=pStop.before & b==3, i) = 2;
        state(temp>=pStop.d_o & d_Inside==3, i) = 2;
        state(temp>=pStop.d_i & d_Outside==3, i) = 2;
        nTurn(state(:,i)==2,i) = nTurn(state(:,i)==2,i-1)+1;                % add to turns counter
        
        % Flies leave the boundary by initiating a curved walk
        state(state(:,i-1)==4 & newTrack, i) = 3;

        % if the fly is at the boundary, then it's in the boundary state
        state(bound & newTrack,i) = 4;
        newTrackState = state(:,i);
        % update the states that are not new tracks
        state(state(:,i)==0,i) = state(state(:,i)==0,i-1);
    end
    
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % set up routines
    
    % cases 1-9 follow the new track routines
    cases = zeros(totFly,15);
    k = 1;
    for s = 1:3
        cases(:,k) = (newTrackState==s) & before;
        cases(:,k+1) = (newTrackState==s) & during & outside;
        cases(:,k+2) = (newTrackState==s) & during & inside;
        k = k+3;
    end
    
    % cases 10-11 follow the new boundary routines
    cases(:,10) = (newTrackState==4) & before;
    cases(:,11) = (newTrackState==4) & during;
    
    % cases 12-14 follow the continue movement routine
    cases(:,12) = ~newTrack & bound;
    cases(:,13) = ~newTrack & state(:,i)<3; % stop and sharp turn have same continue movement routine
    cases(:,14) = ~newTrack & state(:,i)==3;% curved walk has own routine
    
    % case 15 follow the move away from edge routine
    if i>1 
        if any(loc(:,i-1)==1)
            cases(:,15) = newTrack & bound & loc(:,i-1)==1;
        end
    end
    %----------------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    % compute new routines (1-11)
    xyPos = [xPos(:,i),yPos(:,i)];
    % new track routines
    for c = 1:9
        if any(cases(:,c))
            cc = cases(:,c)==1;
            [sampDur,sampSpeed,currAng,sampCurv] = newTrackRoutine...
                (PDFs,turnInProb,currAng,cc,c,xyPos,turningBias,params,dDur,options,nTurn(:,i));
            currDur(cc) = sampDur(cc)-1;
            oldDur(cc) = sampDur(cc);
            currSpd(cc) = sampSpeed(cc);
            currCurv(cc) = sampCurv(cc);
        end
    end
    allRout1 = any(cases(:,1:9),2);
    xPos(allRout1,i+1) = xPos(allRout1,i)+currSpd(allRout1).*cosd(currAng(allRout1));
    yPos(allRout1,i+1) = yPos(allRout1,i)+currSpd(allRout1).*sind(currAng(allRout1));
    
    % new boundary routines
    for c = 10:11
        if any(cases(:,c))
            cc = cases(:,c)==1;
            [sampDur,sampPhi] = boundaryRoutine(PDFs.all,currAng,totFly,c);
            currPhi(cc) = sampPhi(cc);
            currDur(cc) = sampDur(cc)-1;
            oldDur(cc) = sampDur(cc);
        end
    end
    allRout2 = find(cases(:,10) | cases(:,11));
    if ~isempty(allRout2)
        for j = 1:length(allRout2)
            k = allRout2(j);
            theta = currPhi(k);
            Arot = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
            xy = Arot * [xPos(k,i); yPos(k,i)];
            xPos(k,i+1) = xy(1);
            yPos(k,i+1) = xy(2);
        end
    end
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % continue movement routines
    % boundary condition
    if any(cases(:,12))
        c = find(cases(:,12));
        for j = 1:length(c)
            k = c(j);
            theta = currPhi(k);
            Arot = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
            xy = Arot * [xPos(k,i); yPos(k,i)];
            xPos(k,i+1) = xy(1);
            yPos(k,i+1) = xy(2);
        end
        currDur(c) = currDur(c)-1;
    end
    % sharp turn/stop
    if any(cases(:,13))
        c = cases(:,13)==1;
        halfDur = (floor(oldDur./2)+1);
        currAng(c & currDur==halfDur) = currAng(c & currDur==halfDur)+currCurv(c & currDur==halfDur);
        xPos(c,i+1) = xPos(c,i)+currSpd(c).*cosd(currAng(c));
        yPos(c,i+1) = yPos(c,i)+currSpd(c).*sind(currAng(c));
        currDur(c) = currDur(c)-1;
    end
    % curved walk
    if any(cases(:,14))
        c = cases(:,14)==1;
        if curvedRun
            currAng(c) = currAng(c)+currCurv(c);
        end
        if noise
            currAng(c) = currAng(c)+0.5.*randn(1);
        end
        xPos(c,i+1) = xPos(c,i)+currSpd(c).*cosd(currAng(c));
        yPos(c,i+1) = yPos(c,i)+currSpd(c).*sind(currAng(c));
        currDur(c) = currDur(c)-1;
    end
    %----------------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    % move away from edge routine
    if any(cases(:,15))
        c = cases(:,15)==1;
        xPos(c,i+1) = xPos(c,i)./1.02;
        yPos(c,i+1) = yPos(c,i)./1.02;
        
        currAng(c) = myatan(xPos(c,i+1),yPos(c,i+1),'degrees',2)+180+20*(rand(1)-0.5);
        currDur(c) = 0;
        oldDur(c) = 0;
        turn(c) = false;
    end
    %----------------------------------------------------------------------
    
    allAng(:,i+1) = currAng;
    % update x,y, and r position for points that move out of the arena
    % bounds
    xTemp = xPos(:,i+1);
    yTemp = yPos(:,i+1);
    
    rTemp = sqrt(xPos(:,i+1).^2+yPos(:,i+1).^2);
    rTemp(rTemp<1) = 1;
    
    xPos(:,i+1) = xTemp./rTemp;
    yPos(:,i+1) = yTemp./rTemp;
    rPos(:,i+1) = sqrt(xPos(:,i+1).^2+yPos(:,i+1).^2);
    
    newTrack(:) = false;
    progressbar(i/(len-1))
end

firstEntry = ones(totFly,1)*len;
for i = 1:totFly
    tmpFE = find(rPos(i,ceil(len/2):len)<border,1);
    if ~isempty(tmpFE)
        firstEntry(i) = tmpFE+floor(len./2)-1;
    end
    
    stops = find(diff(xPos(i,:))==0 & diff(yPos(i,:))==0);
    
    xyStopsBefore = [xPos(i,stops(stops<firstEntry(i)));yPos(i,stops(stops<firstEntry(i)))];
    xyStopsDuring = [xPos(i,stops(stops>=firstEntry(i)));yPos(i,stops(stops>=firstEntry(i)))];
    
    [uniqueStops,~,ib] = unique(xyStopsBefore', 'rows');
    numoccurences = accumarray(ib, 1);
    uniqueStops = uniqueStops(numoccurences>stopLength,:)';
    synthFlys.stopsBefore{i} = uniqueStops;
    
    [uniqueStops,~,ib] = unique(xyStopsDuring', 'rows');
    numoccurences = accumarray(ib, 1);
    uniqueStops = uniqueStops(numoccurences>stopLength,:)';
    synthFlys.stopsDuring{i} = uniqueStops;
    
end

synthFlys.x = xPos.*4;
synthFlys.y = yPos.*4;
synthFlys.r = rPos.*4;
synthFlys.firstEntry = firstEntry;
end

function [sampDur,sampSpeed,currAng,sampCurv] = newTrackRoutine...
    (PDFs,turnInProb,currAng,caseNdx,i,xyPos,turningBias,params,dDur,options,nTurn)

currPDFs = PDFs.all;
pTurnIn = turnInProb;
totFly = params.totFly;
radPos = sqrt(xyPos(:,1).^2+xyPos(:,2).^2);

sampVals = defineVars(currPDFs.val{i},currPDFs.prob{i},totFly);
%dDur(i) = 3;%diff(unique(currPDFs.val{i}(:,1)));

if size(sampVals,2)==3
    % curved walks
    sampSpeed = sampVals(:,2);
    sampDur = sampVals(:,1)+dDur(i).*rand(totFly,1);
    sampCurv = sampVals(:,3);%+0.05;
    sampCurv(sampCurv.*sampDur>180) = 180./(sampDur(sampCurv.*sampDur>180));
    
    Flys = find(caseNdx);
    [curvDir] = walkTurnBias(xyPos,currAng,sampCurv,sampSpeed,sampDur,...
        radPos(Flys),pTurnIn.all{i},Flys,turningBias(Flys));
    sampCurv(caseNdx) = curvDir;
else
    % sharp turn and stops
    sampDur = sampVals(:,1)+dDur(i).*rand(totFly,1);
    sampCurv = sampVals(:,2);
    [curvDir] = turnBias(xyPos,currAng,sampCurv,pTurnIn.all{i},radPos,turningBias);
    sampCurv = curvDir;
end
if i<4
    % stops
    sampSpeed = zeros(totFly,1);
elseif i<7
    % sharp turns
    sampSpeed = 0.*ones(totFly,1);%0.02.*ones(totFly,1);%0.03.*ones(totFly,1);
end
sampDur = floor(sampDur)+1;
cf = 40;% 40 mm for radius
% if i >7 && options.borderChoice
%     cf(nTurn<3) = 40;
%     cf(nTurn>=3) = 40;
%     %cf = 35;
% end
sampSpeed = sampSpeed./cf;

% sharp turn
if i > 6
    currAng(caseNdx) = currAng(caseNdx)+sampCurv(caseNdx);
end

end

function [synthVar] = defineVars(val,prob,totFly)

x = discretesample2(prob,totFly);
synthVar = val(x,:);

end

function [turnDir] = turnBias(xyPos,currAng,sampYaw,pTurnIn,radPos,turningBias)
x = xyPos(:,1); y = xyPos(:,2);
nFlys = size(x,1);
%u = [cosd(currAng),sind(currAng),zeros(nFlys,1)];
v1 = [cosd(currAng+sampYaw),sind(currAng+sampYaw),zeros(nFlys,1)];
v2 = [cosd(currAng-sampYaw),sind(currAng-sampYaw),zeros(nFlys,1)];
w = [-x,-y,zeros(nFlys,1)];

% calculate the angles between the center facing vector and the
% direction the fly is moving in before the turn and after the turn
a = cross(v1,w);
b = dot(v1',w')';
c = sqrt(sum(a.^2,2));
turnLeft = atan2(c,b).*180./pi;
a = cross(v2,w);
b = dot(v2',w')';
c = sqrt(sum(a.^2,2));
turnRight = atan2(c,b).*180./pi;

% if any(turnRight<0) || any(turnLeft<0)
%     a
% end

turnDir = sampYaw;
turnDir(turnRight<turnLeft) = -sampYaw(turnRight<turnLeft);

[sigTemp] = computeDir(radPos,pTurnIn,turningBias);
turnDir(sigTemp) = turnDir(sigTemp).*-1;

end

function [curvDir] = walkTurnBias(xyPos,currAng,sampCurv,sampSpeed,sampDur,radPos,pTurnIn,Flys,turningBias)
x = xyPos(:,1); y = xyPos(:,2);
nFlys = length(Flys);
xEnd = zeros(nFlys,2);yEnd = zeros(nFlys,2);
for f = 1:nFlys
    j = Flys(f);
    xEnd(f,1) = x(j)+sum(sampSpeed(j).*cosd(currAng(j)+sampCurv(j).*[1:1:sampDur(j)]));
    yEnd(f,1) = y(j)+sum(sampSpeed(j).*cosd(currAng(j)+sampCurv(j).*[1:1:sampDur(j)]));
    xEnd(f,2) = x(j)+sum(sampSpeed(j).*cosd(currAng(j)-sampCurv(j).*[1:1:sampDur(j)]));
    yEnd(f,2) = y(j)+sum(sampSpeed(j).*cosd(currAng(j)-sampCurv(j).*[1:1:sampDur(j)]));
end

turnLeft = sqrt((xEnd(:,1)).^2+(yEnd(:,1)).^2);
turnRight = sqrt((xEnd(:,2)).^2+(yEnd(:,2)).^2);

curvDir = sampCurv(Flys);
curvDir(turnRight<turnLeft) = -sampCurv(turnRight<turnLeft);

[sigTemp] = computeDir(radPos,pTurnIn,turningBias);

curvDir(sigTemp) = curvDir(sigTemp).*-1;

end

function [dir] = computeDir(rPos,pTurnIn,turningBias)

sBLim = (0:0.2:size(pTurnIn,2)*0.2)./4;sBLim(end) = sBLim(end)+0.1;
if sum((rPos<sBLim(end) & rPos>=sBLim(1)))>0
    pSamp = rand(length(rPos),1);
    dir = false(size(pSamp));
    if all(pTurnIn == pTurnIn(1,1))
        dir = pSamp>=pTurnIn(1,1);
    else
        for j = 1:size(pTurnIn,2)
            bordLoc = (rPos<sBLim(j+1) & rPos>=sBLim(j));
            SampChoice = pSamp>pTurnIn(1,j);
            SampChoice2 = pSamp>pTurnIn(2,j);
            SampChoice3 = pSamp>0.5;
            
            assert(length(bordLoc)==length(SampChoice),'err1')
            assert(length(bordLoc)==length(turningBias),'err2')
            
            dir(bordLoc & SampChoice & turningBias==1) = true;
            dir(bordLoc & SampChoice2 & turningBias==2) = true;
            dir(bordLoc & SampChoice3 & turningBias==0) = true;
        end
    end
end
end

function [sampDur,sampPhi] = boundaryRoutine(stopWalkPDFs,currAng,totFly,i)

synthVar = defineVars(stopWalkPDFs.val{i},stopWalkPDFs.prob{i},totFly);
%sampDur = floor(synthVar(:,1))+1;
sampDur = floor(synthVar(:,1))+10;
sampPhi = synthVar(:,2);
sampPhi = sampPhi./sampDur;
sampPhi = sampPhi.*sign(currAng);

end


