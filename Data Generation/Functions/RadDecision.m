function [durProbDetrend,durProb,turnInDuringProb,ratOutIn,ratInOut,turnAmount] = RadDecision(matFile,matFileBefore,opts)
% This function calculates turn bias and border choice
%
% Inputs:
%    matFile: file name for crossing file following first entry
%    matFileBefore: file name for crossing file before first entry
%    opts: structure with field border

% Outputs:
%    durProbDetrend: 2xn matrix of detrended turn in probability
%       row 1: following first entry 
%       row 2: before first entry
%    durProb: 2xn matrix of undetrended turn in probability
%       row 1: following first entry 
%       row 2: before first entry
%    turnInDuringProb: 2xn matrix of probability of turning in
%       row 1: first 2 turns
%       row 2: later turns
%    ratOutIn: mx1 matrix of percent of tracks that display at least m
%    turns after entering the light ring following first entry
%    ratInOut: mx1 matrix of percent of tracks that display at least m
%    turns after leaving the light ring following first entry
%    turnAmount: Structure with fields
%       turnAmount.first2TurnsOutIn: amount turned for first two turns
%       after entering light ring
%       turnAmount.first2TurnsInOut: amount turned for first two turns
%       after leaving light ring
%       turnAmount.laterTurnsOutIn: amount turned for later turns after
%       entering light ring
%       turnAmount.laterTurnsInOut: amount turned for later turns after
%       leaving light ring

warning('off','curvefit:fit:noStartPoint');
close all
% This was the empirically measured intensity (for plotting purposes only)
Intensity = [26030, 26020, 26090, 25630, 25800, 25440, 25500, 24920, ...    % Measured intensity in nW
25030, 24780, 25100, 16700, 3610, 134, 31.5, 20];
I = [Intensity./max(Intensity) zeros(1,24)];
border = opts.border;

figure;
% Generate border choice and turn bias for first 2 turns before first entry
fName = matFileBefore;
[before_totDec,turnInBeforeProb,~,~,~,~,~] = generateDistributions(fName,I,border,[1:2],1);
% Generate border choice and turn bias for first 2 turns after first entry
fName = matFile;
[during_totDec,turnInDuringProb,x,ratOutInDur,ratInOutDur,~,turnTmp] = generateDistributions(fName,I,border,[1:2],3);
set(gcf,'Position',[9 49 1014 918])
turnAmount.first2TurnsOutIn = turnTmp.OutIn;
turnAmount.first2TurnsInOut = turnTmp.InOut;

figure;
% Generate border choice and turn bias for later turns before first entry
fName = matFileBefore;
[before_totDec2,turnInBeforeProb2,~,~,~,~,~] = generateDistributions(fName,I,border,[3:100],1);
% Generate border choice and turn bias for later turns after first entry
fName = matFile;
[during_totDec2,turnInDuringProb2,~,ratOutInDur2,ratInOutDur2,~,turnTmp] = generateDistributions(fName,I,border,[3:100],3);
set(gcf,'Position',[9 49 1014 918])
turnAmount.laterTurnsOutIn = turnTmp.OutIn;
turnAmount.laterTurnsInOut = turnTmp.InOut;

% set unknown turn bias to 50/50
% turnInBeforeProb(isnan(turnInBeforeProb)) = 0.5;
% turnInDuringProb(isnan(turnInDuringProb)) = 0.5;
% turnInBeforeProb2(isnan(turnInBeforeProb2)) = 0.5;
% turnInDuringProb2(isnan(turnInDuringProb2)) = 0.5;

% percent of tracks with at least n turns 
ratOutIn = [ratOutInDur;ratOutInDur2];
ratInOut = [ratInOutDur;ratInOutDur2];

% get the detrend function from before first entry
ratio = sum(during_totDec)./sum(before_totDec);
y = before_totDec'.*ratio;
x1 = x(1:14);
y1 = y(1:14);
expEqn = 'u.*exp(u.*(x-a))';
fLeft = fit(x1',y1,expEqn);
% model In to Out
x2 = x(14:29);
y2 = y(14:29);
expEqn = 'u.*exp(-u.*(x-a))';
fRight = fit(x2',y2,expEqn);
trendLeft = fLeft.u.*exp((fLeft.u).*(x1-fLeft.a));
trendRight = fRight.u.*exp(-(fRight.u).*(x(15:end)-fRight.a));
trend1 = [trendLeft trendRight];

% subtract out the detrend function
detrend_During_dec = during_totDec-trend1;
detrend_During_dec(detrend_During_dec<0) = 0;

% later turns do not need the detrend function
trend1 = [zeros(size(during_totDec2))];
detrend_During_dec2 = during_totDec2-trend1;
detrend_During_dec2(detrend_During_dec2<0) = 0;

% consolidate to single array
turnInDuringProb = [turnInDuringProb;turnInDuringProb2];
durProbDetrend = [detrend_During_dec;detrend_During_dec2];
durProb = [during_totDec;during_totDec2];

durProbDetrend(isnan(durProbDetrend)) = 0;
%turnInDuringProb(isnan(turnInDuringProb)) = 0.5;

% plot turn bias for the first 2 turns
figure;subplot(2,1,1);plot(x(1:2:end),turnInBeforeProb);
xlim([0 1]);title(' Turn in Prob Before');ylim([0 1])
ylabel('Probability');xlabel('r Distance');suptitle('First 2 turns')
subplot(2,1,2);plot(x(1:2:end),turnInDuringProb(1,:));title(' Turn in Prob During')
hold on;plot([1.3/4 1.3/4],[0 1],'r');xlim([0 1]);
ylabel('Probability');xlabel('r Distance');suptitle('First 2 turns')

% plot turn bias for consecutive turns
figure;subplot(2,1,1);plot(x(1:2:end),turnInBeforeProb2);
xlim([0 1]);title(' Turn in Prob Before')
ylabel('Probability');xlabel('r Distance');suptitle('Later turns')
subplot(2,1,2);plot(x(1:2:end),turnInDuringProb2);title(' Turn in Prob During')
hold on;plot([1.3/4 1.3/4],[0 1],'r');xlim([0 1]);
ylabel('Probability');xlabel('r Distance');suptitle('Later turns')

% plot the probability of turning (normalize for inside) and for outside
% separately
figure;plot([1:12]./40,during_totDec2(1:12)./sum(during_totDec2(1:12)),...
    '--','Color',[0.5 0.5 0.5]);hold on;
plot([13:40]./40,during_totDec2(13:40)./sum(during_totDec2(13:40)),...
    '--','Color',[0.5 0.5 0.5]);
plot([1:12]./40,during_totDec(1:12)./sum(during_totDec(1:12)),'k');hold on;
plot([13:40]./40,during_totDec(13:40)./sum(during_totDec(13:40)),'k');
ylim([0 0.4])
legend({'later turns','later turns','first 2 turns','first 2 turns'})

end

function [totDec,turnInProb,x,ratOutIn,ratInOut,ProbNumTurn,turnAmount] = generateDistributions(fName,I,border,n,subplt)
load(fName,'empFlys','curvPks','tmpOutIn','tmpInOut','tmpTracksOutIn2','tmpTracksInOut2')
% radius of arena is 4 cm
x = (0.1:0.1:4.1)/4;
empFlys.rH = empFlys.rH./4;
border = border*4;

% out to in
tmp = tmpOutIn;
tmpTracks2 = tmpTracksOutIn2;
[outInLoc,outInTime,angRelIn,ratOutIn,ProbNumTurn,~,turnAmount.OutIn] = genDist(empFlys,curvPks,tmp,tmpTracks2,I,x,n,subplt);

% in to out
tmp = tmpInOut;
tmpTracks2 = tmpTracksInOut2;
[inOutLoc,inOutTime,angRelIn2,ratInOut,~,~,turnAmount.InOut] = genDist(empFlys,curvPks,tmp,tmpTracks2,I,x,n,subplt);

% plot marginal location to turn pmf in 3d
subplot(2,2,subplt);
dists1 = histcounts(inOutLoc,(border:0.1:4.1)./4);
dists1 = dists1./sum(dists1);
plot3(310*ones(1,length(dists1)),(border:0.1:4)./4,dists1,'b');hold on
dists2 = histcounts(outInLoc,(0:0.1:border)./4);
dists2 = dists2./sum(dists2);
plot3(310*ones(1,length(dists2)),(0:0.1:border-0.1)./4,dists2,'k');

% plot light boundary in 3d
plot3(310*ones(1,length(I)),(0:0.1:3.9)./4,I./3,'r');

% plot marginal time to turn pmf in 3d
dists3 = histcounts(inOutTime,(0:10:320));
dists3 = dists3./sum(dists3);
plot3((0:10:310),ones(1,length(dists3)),dists3,'b');
dists4 = histcounts(outInTime,(0:10:320));
dists4 = dists4./sum(dists4);
plot3((0:10:310),ones(1,length(dists4)),dists4,'k');
view(3)
xlabel('Time');ylabel('Radial Distance');zlabel('Probability')

totDec = histcounts(inOutLoc,(0:0.1:4.1)./4)+histcounts(outInLoc,(0:0.1:4.1)./4);
% calculate turns where the fly turns inwards
allDec = [outInLoc,inOutLoc];
angRelInAll = [angRelIn,angRelIn2];
turnInLoc = angRelInAll(1,:)>angRelInAll(2,:);

tmp = [histcounts(allDec(turnInLoc==1),(0:0.2:4.2)./4);histcounts(allDec,(0:0.2:4.2)./4)];
turnInProb = tmp(1,:)./tmp(2,:);turnInProb(tmp(2,:)<5) = nan;

% plot border choice
subplot(2,2,subplt+1);hold on
plot(x,histcounts(allDec(turnInLoc==1),(0:0.1:4.1)./4));hold on;plot(x,histcounts(allDec,(0:0.1:4.1)./4))
legend({'Turns Inwards','All Turns'})
xlabel('r Distance');ylabel('# of Instances')

end


function [outInLoc,outInTime,angRelIn,ratOutIn,ProbNumTurn,NTurns,turnAmount] = genDist(empFlys,curvPks,tmp,tmpTracks2,I,x,n,subplt)

tmp2 = [];angRelIn = [];
% loop through each crossing track
for j = 1:size(tmp,1)
    % get the fly and track number of the track
    currFly = tmp(j,3)+1;
    currTrack = tmpTracks2(j,:);
    currTrack = currTrack(~isnan(currTrack));
    currTrack(currTrack==0) = [];
    if ~isempty(currTrack)
        % find all sharp turns that occur within range of the track
        currNdx = find(currTrack(1)<curvPks.ndx{1,currFly} & curvPks.ndx{1,currFly}<currTrack(end));
        if ~isempty(currNdx)
            % when the sharp turn happens
            SharpTurnAbsNdx = curvPks.ndx{1,currFly}(currNdx);
            % when the sharp turn happens relative to the start of the
            % track
            SharpTurnRelNdx = SharpTurnAbsNdx-currTrack(1);
            % head position at sharp turns
            SharpTurnRLoc = empFlys.rH(currFly,SharpTurnAbsNdx);
            % orientation before and after sharp turn
            angBef = curvPks.dirRelCenterBef{1,currFly}(currNdx);
            angAft = curvPks.dirRelCenterAft{1,currFly}(currNdx);
            % total curvature of sharp turn
            curvAmount = curvPks.tot{1,currFly}(currNdx);
            
            curvDistance = nan(1,numel(currNdx));
            curvAmount2 = nan(1,numel(currNdx));
            a = nan(1,numel(currNdx));
            for k = 1:numel(currNdx)
                tmpNdx = curvPks.all{1,currFly}{1,currNdx(k)}(:,2);
                % curve distance and curve amount
                curvDistance(k) = sum(sqrt(empFlys.thrust(currFly,tmpNdx).^2+empFlys.slip(currFly,tmpNdx).^2));
                curvAmount2(k) = sum(curvPks.all{1,currFly}{1,currNdx(k)}(:,1));
                
                xTmp = empFlys.x(currFly,tmpNdx);xTmp = [xTmp xTmp(1)];
                yTmp = empFlys.y(currFly,tmpNdx);yTmp = [yTmp yTmp(1)];
                % area covered by the sharp turn (loopyness of curve)
                a(k) = polyarea(xTmp',yTmp');
            end
            
            % consolidate parameters into arrays
            angRelIn = [angRelIn,[angBef;angAft]];
            tmp2 = [tmp2, [SharpTurnAbsNdx;SharpTurnRelNdx;SharpTurnRLoc;(1:1:length(SharpTurnAbsNdx));curvAmount;curvDistance;a]];
        end
    end
end

% loop through each turn index (i.e first turn, second turn, etc since
% crossing)
for j = 1:length(n)
    firstTurn(j,:) = tmp2(4,:)==n(j);
    nTurn(j,:) = tmp2(4,:)==n(j);
end
ratOutIn = sum(firstTurn,2)./size(tmpTracks2,1);
firstTurn = any(firstTurn,1);
nTurn = any(nTurn,1);

% assign interested parameters
outInLoc = tmp2(3,firstTurn);
outInTime = tmp2(2,firstTurn);
angRelIn = angRelIn(:,firstTurn);
turnAmount = tmp2([5,6,7],firstTurn);
ProbNumTurn{1} = histcounts(tmp2(4,:),[1:1:max(tmp2(4,:))])./size(tmpTracks2,1);
NTurns.outIn = size(tmpTracks2,1);

% plot 85th percentile contour for joint pmf of turn location and time
xT = tmp2(2,nTurn);
yT = tmp2(3,nTurn);
if ~isempty(xT)
    [pdf,X1,X2]=akde([xT',yT']);
    pdf = pdf./sum(pdf);

    m=sort(pdf,'descend');
    v=max(find(cumsum(m)<0.85));
    pVal=m(v);

    pdfAll=reshape(pdf,size(X1));

    X = X1(1,:);
    Y = X2(:,1)';

    subplot(2,2,subplt);hold on
    c = parula(5);
    for j = 1:length(x)-1
        fill([0 0 300 300],[x(j) x(j+1) x(j+1) x(j)],'r','EdgeColor','none','FaceAlpha',I(j)./2)
    end
    contour(X,Y,pdfAll,[pVal pVal],'LineWidth',2,'Color','k')
        plot(tmp2(2,firstTurn),tmp2(3,firstTurn),'.','Color','k','markersize',10);
    xlim([0 310]);ylim([0 1])
end

end



