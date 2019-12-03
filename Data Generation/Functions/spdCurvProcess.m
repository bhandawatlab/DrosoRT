function [rad,speed,curv,allScenarios,phi,pStop,empFlys] = spdCurvProcess(data,nFlies,fs,spdThresh,boundThresh)
% This function processes the separates out the fly's trajectory into one
% of 8 scenarios:
%       before stop; during Outside stop; during Inside stop; before run; 
%       during Outside run; during Inside run; before boundary; during boundary;
% This is to be used downstream in calculating joint pmfs for stop, curved
% walks, sharp turns, and boundary states
%
% Inputs:
%    data: Cell array of data (1 x n Fly)
%       for each cell, in descending rows:
%           x;y;r;in;during;thrust;slip;yaw;curvFly;xH;yH;rH
%    nFlies: number of flies
%    fs: sampling frequency
%    spdThresh: speed threshold for stopping
%    boundThresh: radial position threshold for boundary condition

% Outputs:
%    rad: radial position (m time points x n flies)
%    speed: speed (m time points x n flies)
%    rad: curvature (m time points x n flies)
%    allScenarios: 3d matrix of scenarios (m time points x n flies x 
%       scenario number)
%    phi: radial angle relative to positive x-axis (m time points x n flies)
%    pStop: structure showing probability of entering a stop
%    empFlys: structure that contrains information about all flies in the
%       genotype

% threshold for plotting purposes to visualize where long stops occur
stopLength = 90;

% obtain all possible speed and curvatures (distribution)
npts = length(data{1})+1;
speed = zeros(npts,nFlies);yaw = zeros(npts,nFlies);%boundary = false(npts,nFlies);
phi = zeros(npts,nFlies);curv = zeros(npts,nFlies);rad = zeros(npts,nFlies);

allScenarios = false(npts,nFlies,8);
% loop through each fly
for i = 1:nFlies
    % get speed, yaw, curvature, and radial position
    spdTemp = sqrt(data{i}(6,:).^2+data{i}(7,:).^2);
    speed(:,i) = [spdTemp 0]'./fs;
    yaw(:,i) = [data{i}(8,:) 0]';
    curv(:,i) = [data{i}(9,:) 0]';
    %rad(:,i) = [data{i}(3,:) 0]';                                          % use radial position defined by body
    rad(:,i) = [data{i}(12,:) 0]';                                          % use radial position defined by head
    
    % index when boundary condition applies
    boundary = (data{i}(3,:)>boundThresh)';
    
    % fill in empFlys structure with data from Data
    empFlys.x(i,:) = data{i}(1,:);
    empFlys.y(i,:) = data{i}(2,:);
    empFlys.r(i,:) = data{i}(3,:);
    empFlys.thrust(i,:) = data{i}(6,:);
    empFlys.slip(i,:) = data{i}(7,:);
    
    empFlys.xH(i,:) = data{i}(10,:);
    empFlys.yH(i,:) = data{i}(11,:);
    empFlys.rH(i,:) = data{i}(12,:);
    
    % define when flies are inside or outside
    in = data{i}(4,:)';
    out = ~in;
    
    % define before/after light on
    dur = data{i}(5,:)';
    bef = ~dur;
    
    % define first entry
    fe = find(dur & in,1);
    empFlys.firstEntry(i,1) = fe;
    FETemp = fe;
    
    % index when flies are stopped or not stopped
    stop = (spdTemp<(spdThresh))';%.*7.7345e-04
    run = ~stop;
    
    % populate the allScenarios matrix with each of the 8 scenarios
    allScenarios(:,i,1) = [bef & ~boundary & stop; false];       % before stop;
    allScenarios(:,i,2) = [dur & out & ~boundary & stop; false]; % during Outside stop;
    allScenarios(:,i,3) = [dur & in & ~boundary & stop; false];  % during Inside stop;
    
    allScenarios(:,i,4) = [bef & ~boundary & run; false];        % before run;
    allScenarios(:,i,5) = [dur & out & ~boundary & run; false];  % during Outside run;
    allScenarios(:,i,6) = [dur & in & ~boundary & run; false];   % during Inside run;
    allScenarios(:,i,7) = [bef & boundary; false];               % before boundary;
    allScenarios(:,i,8) = [dur & boundary; false];               % during boundary;
    
    % angle from the positive x-axis
    phi(:,i) = [myatan(data{i}(1,:),data{i}(2,:),'degrees',2),0]';
    
    % this part is just to assign long stops as a field in empFlys. Use
    % this field for visualization purposes only
    startStop = find(diff([0; stop])==1);
    endStop = find(diff([stop;0])==-1);
    numoccurences = endStop-startStop;
    startStop = startStop(numoccurences>stopLength)+30;
    xyStopsBefore = [empFlys.x(i,startStop(startStop<FETemp));empFlys.y(i,startStop(startStop<FETemp))];
    xyStopsDuring = [empFlys.x(i,startStop(startStop>=FETemp));empFlys.y(i,startStop(startStop>=FETemp))];
    empFlys.stopsBefore{i} = xyStopsBefore;
    empFlys.stopsDuring{i} = xyStopsDuring;
end

% define probability of transitioning into a stop for before, during
% outside, and during inside
pStop.before = sum(sum(allScenarios(:,:,1)))./(sum(sum(allScenarios(:,:,1)))+sum(sum(allScenarios(:,:,4))));
pStop.d_o = sum(sum(allScenarios(:,:,2)))./(sum(sum(allScenarios(:,:,2)))+sum(sum(allScenarios(:,:,5))));
pStop.d_i = sum(sum(allScenarios(:,:,3)))./(sum(sum(allScenarios(:,:,3)))+sum(sum(allScenarios(:,:,6))));
end
