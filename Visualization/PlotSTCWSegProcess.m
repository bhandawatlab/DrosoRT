function [] = PlotSTCWSegProcess(GlobMinX)
% This function plots an example long segment of a fly track to show what
% parts were delineated into sharp turns and curved walks
%
% Inputs:
%    GlobMinX: best fit parameters for delineating between sharp turns and
%       curved walks

% load curve walk, sharp turn, and data matrices
load([pwd '\Data Full\DataCons\Orco Retinal_Nov13.mat'],'curvPks','curvWalks','Data');
% load parameter fits
close all

% specify parameters
fly = 2;
fs = 30;
nST = length(curvPks.all{fly});
nCW = length(curvPks.all{fly});
len = size(Data.curv,2);

% create type vectors to mark locations of sharp turn and curved walks
type = zeros(1,len);type2 = zeros(1,len);
for i = 1:nST
    type(curvPks.all{fly}{i}(:,2))=1;
    type2(curvPks.all{fly}{i}(:,2))=1;
end
for i = 1:nCW
    type(curvWalks.all{fly}{i}(:,2))=2;
    type2(curvWalks.all{fly}{i}(:,2))=1;
end

% find when there are stops or boundary conditions and separate the track
% into trajectories of only sharp turn and curved walk
[startNdx,endNdx,ndx] = startEndSeq(type2);
startNdx(ndx==0) = [];
endNdx(ndx==0) = [];

% choose the longest trajectory
trackLen = endNdx-startNdx;
track = [startNdx(trackLen==max(trackLen)):endNdx(trackLen==max(trackLen))];
trackType = type(track);

% find start and end indexes of sharp turns within the chosen trajectory
[startNdx,endNdx,ndx] = startEndSeq(trackType);
startNdx(ndx==2) = [];
endNdx(ndx==2) = [];

% calculate the thresholds in degrees/frame
cThresh = pi/(GlobMinX(4)*400).*180./pi;
dcThresh = pi/(GlobMinX(5)*400).*180./pi;

% create new variables for the curvature and change in curvature tracks
curvTrack = abs(Data.curv(fly,track)).*180./pi;
% dcurv = abs(diff(Data.curv(fly,:)));
% dcurvTrack = dcurv(track).*180./pi;
dcurvTrack = [abs(diff(curvTrack)) 0];

% specify the weighted values sequence of binarized curvature and change in
% curvature
weight_Curv = GlobMinX(1).*(curvTrack>cThresh);
weight_dCurv = GlobMinX(2).*(dcurvTrack>dcThresh);
TFlag = weight_Curv+weight_dCurv;
sTFlag = smooth(TFlag,ceil(GlobMinX(3)*16));
turnThresh = (1+0.6)*mean(sTFlag);
walkThresh = (1-0.6)*mean(sTFlag);

% specify x in seconds
x = (track-track(1)+1)./fs;

% plot the curvature of the trajectory with binarized curvature in red.
% In addition, plot the curvature threshold parameter
figure;subplot(5,1,1)
plot(x,curvTrack,'k','linewidth',1);hold on
plot([0 max(x)],[cThresh cThresh],'-b','linewidth',1)
xlim([0 x(end)])
xlabel('time (s)');ylabel('Curv (deg/frame)')

% plot the change in curvature of the trajectory with binarized 
% In addition, plot the change in curvature threshold parameter
subplot(5,1,2)
plot(x,dcurvTrack,'k','linewidth',1);hold on
plot([0 max(x)],[dcThresh dcThresh],'-g','linewidth',1)
xlim([0 x(end)])
xlabel('time (s)');ylabel('dCurv (deg/(frame^2))')

subplot(5,1,3)
stairs(x,TFlag,'-g','linewidth',1);hold on
stairs(x,weight_Curv,'-b','linewidth',1)
xlim([0 x(end)])
xlabel('time (s)');ylabel('Weighted sum')

subplot(5,1,4);hold on
stairs(x,sTFlag,'Color',[0.5 0.5 0.5],'Linewidth',1);
plot([0 max(x)],[turnThresh turnThresh],'Color','r','Linewidth',1)
plot([0 max(x)],[1 1]*(walkThresh+turnThresh)./2,'Color',[1 0.4 0.6],'Linewidth',1)
plot([0 max(x)],[walkThresh walkThresh],'Color','k','Linewidth',1)
xlim([0 x(end)])
xlabel('time (s)');ylabel('Smoothed weighted sum')

subplot(5,1,5);hold on
plot(x,curvTrack,'k','linewidth',1);hold on
for i = 1:length(startNdx)
    xtmp = (startNdx(i)-1:endNdx(i)+1)./fs;
    plot(xtmp,curvTrack(startNdx(i)-1:endNdx(i)+1),'r','linewidth',1);
end
xlim([0 x(end)])
xlabel('time (s)');ylabel('Curv (deg/frame)')
suptitle('Sample ST/CW Trajectory')
set(gcf,'Position',[1.6667, 41.6667, 638.6667, 599.3333])

% plot the xy positions of the trajectory with sharp turn regions as
% designated by red. In addition, plot the arena boundary in gray and the
% rest of the fly track in green
figure;
plotCircle([0 0],4,100,[0.5 0.5 0.5]);hold on
plot(Data.x(fly,:),Data.y(fly,:),'g','linewidth',2)
plot(Data.x(fly,track),Data.y(fly,track),'k','linewidth',2)
for i = 1:length(startNdx)
    plot(Data.x(fly,track(startNdx(i):endNdx(i))),Data.y(fly,track(startNdx(i):endNdx(i))),'r','linewidth',2);
end
axis equal;xlim([-10 5])
text(-9,3,['Start Pos: (' num2str(Data.x(fly,track(1))) ', ' ...
    num2str(Data.y(fly,track(1))) ')'],'FontSize',8)
text(-9,2,[' End Pos: (' num2str(Data.x(fly,track(end))) ', ' ...
    num2str(Data.y(fly,track(end))) ')'],'FontSize',8)
suptitle('Sample ST/CW Trajectory')
set(gcf,'Position',[1.6667, 41.6667, 638.6667, 599.3333])


end


