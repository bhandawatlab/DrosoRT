function [] = compareDistributions(curvPksKatsov,curvWalksKatsov)
% Compare the kinematic distributions between the Katsov dataset L and our
% dataset
%
% Inputs:
%    curvPksKatsov: structure with information of where sharp turns are in
%       the Katsov tracks
%    curvWalksKatsov: structure with information of where curved walks are 
%       in the Katsov tracks

file1 = [pwd '\Data Full\DataCons\Orco Control_Nov13.mat'];
file2 = [pwd '\Data Full\DataCons\Orco Retinal_Nov13.mat'];
[ST_Cont,CW_Cont] = genOrcoData(file1);
[ST_Ret,CW_Ret] = genOrcoData(file2);

fileKatsov = 'Data Full\KatsovData\dataset_L.mat';
[ST_Katsov,CW_Katsov] = genKatsovData(fileKatsov,curvPksKatsov,curvWalksKatsov);

figure;
type = {'speed','yaw'};
lim = [0 35;-20 20];
edges{1} = [0.3:0.3:35];edges{2} = [-180:0.5:180];
s(1:2) = {CW_Cont};s(3:4) = {CW_Ret};s(5:6) = {CW_Katsov};
s2(1:2) = {ST_Cont};s2(3:4) = {ST_Ret};s2(5:6) = {ST_Katsov};
titles = {'Orco Control Speed','Orco Control Yaw','Orco Retinal Speed',...
    'Orco Retinal Yaw','Katsov Speed','Katsov Yaw'};
plotType = 'count';
for i = 1:6
    n = 2-mod(i,2);
    plotHist(s{i},i,type{n},edges{n},lim(n,:),plotType);hold on;
    plotHist(s2{i},i,type{n},edges{n},lim(n,:),plotType)
    title(titles{i})
end
set(gcf,'Position',[2 42 638 600])
end

function [] = plotHist(s,i,type,edges,lim,plotType)
subplot(3,2,i);
histogram(s.(type)(s.speed>0.3),edges,'Normalization',plotType);
xlim(lim)
end

function [ST,CW] = genOrcoData(file)

load(file,'empFlys','curvPks','curvWalks','Data');
yawAll = Data.yaw;
spdAll = sqrt(Data.thrust.^2+Data.slip.^2);

border = 1/4;% border is arbitrary in this case
flies = length(empFlys.firstEntry);
Data = convertSynthetic(empFlys,flies,border);
nCt = length(Data);

err = zeros(size(spdAll,1),nCt);
l = size(Data{1},2);
for i = 1:size(spdAll,1)
    for j = 1:nCt
        err(i,j) = sum(abs(spdAll(i,1:l)-Data{j}(6,:).*10))./l;
    end
end
globErr = min(err');
num2del = i-j;
[~,ndx] = sort(globErr,'descend');

yawAll(ndx(1:num2del),:) = [];
yawAll = yawAll(:,1:l);

speedST = [];speedCW = [];
yawST = [];yawCW = [];
for i = 1:nCt
    speed = sqrt((Data{i}(6,:).*10).^2+(Data{i}(7,:).*10).^2);
    yaw = yawAll(i,:);
    
    tmpSTNdx = false(size(yaw));tmpCWNdx = false(size(yaw));
    if ~isempty(curvPks.all{1,i})
        tmpST = cell2mat(curvPks.all{1,i}');
        tmpCW = cell2mat(curvWalks.all{1,i}');
        tmpSTNdx(tmpST(:,2)) = true;
        tmpCWNdx(tmpCW(:,2)) = true;
    end
    tmpSTNdx = tmpSTNdx(1:numel(yaw));
    tmpCWNdx = tmpCWNdx(1:numel(yaw));
%    tmpCWNdx = ~tmpSTNdx;
    
    speedST = [speedST speed(tmpSTNdx)];
    speedCW = [speedCW speed(tmpCWNdx)];
    
    yawST = [yawST yaw(tmpSTNdx)];
    yawCW = [yawCW yaw(tmpCWNdx)];
end

ST.speed = speedST;
ST.yaw = yawST;

CW.speed = speedCW;
CW.yaw = yawCW;

%figure;histogram(speedST);hold on;histogram(speedCW);plot([0.3 0.3],[0 6e4],'k')


end

function [ST,CW] = genKatsovData(file,curvPksKatsov,curvWalksKatsov)
load(file,'UI','VR','VRTH','VS','VT')
nCt = max(UI);
fs = 30;

speedST = [];speedCW = [];
yawST = [];yawCW = [];
for i = 1:nCt
    currTrack = find(UI==i);
    
    thrust = VT(currTrack)';
    slip = VS(currTrack)';
    speed = sqrt(thrust.^2+slip.^2);
    yaw = VR(currTrack)'./fs;
    
    tmpSTNdx = false(size(yaw));tmpCWNdx = false(size(yaw));
    if ~isempty(curvPksKatsov{1,i}.all)
        tmpST = cell2mat(curvPksKatsov{1,i}.all');
        tmpCW = cell2mat(curvWalksKatsov{1,i}.all');
        tmpSTNdx(tmpST(:,2)) = true;
        tmpCWNdx(tmpCW(:,2)) = true;
    end
    %tmpCWNdx = ~tmpSTNdx;
    
    speedST = [speedST speed(tmpSTNdx)];
    speedCW = [speedCW speed(tmpCWNdx)];
    
    yawST = [yawST yaw(tmpSTNdx)];
    yawCW = [yawCW yaw(tmpCWNdx)];
    
end

ST.speed = speedST.*10;
ST.yaw = yawST;
CW.speed = speedCW.*10;
CW.yaw = yawCW;
end
