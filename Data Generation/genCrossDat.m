function [] = genCrossDat(opts)
% This function gets crossing tracks and calculates turn bias and border
% choice parameters
%
% Inputs:
%    opts: Type of points we want to return (3 options)
%       opts.genFold = folder with the general files of fly movement 
%           information
%       opts.STCWFold = folder with the sharp turn, curved walk, stop, and
%           boundary conditions information
%       opts.crossFold = folder to save the fly crossing light ring
%           information
%       opts.TBFold = folder to save the fly turn bias information

% get data folders
DataFold = opts.genFold;
STCWFold = opts.STCWFold;
newFold = opts.crossFold;
% multiply by 4 because the arena is 4 mm in radius and all tracks are
% currently in normalized form (max = 1)
border = opts.border.*4;

% get files in the folders
s = dir([DataFold '/*.mat']);
s2 = dir([STCWFold '/*.mat']);
filelist = {s.name};
filelistSTCW = {s2.name};
nGen = length(filelist);

% Loop through each genotype
fnameBefore = cell(1,nGen);fnameAfter = cell(1,nGen);
for i = 1:nGen
    % clear unneccesary variables
    clearvars -except filelist s i filelistSTCW s2 DataFold STCWFold ...
        newFold border opts fnameBefore fnameAfter
    % load data
    load([DataFold '/' filelist{i}],'Data','Arena','fs')
    fliesOld = length(Data.lightOn);
    % convert data to a different format
    [Data,badFlies] = convertFormat(Data,Arena,border,fliesOld);
    flies = length(Data);
    
    % load the four movement states
    load([STCWFold '/' filelistSTCW{i}],'curvPks','curvWalks','stopCond','boundCond')
    
    % remove any flies that did not enter the center light zone in the
    % during case
    [curvPks] = removeBadData(curvPks,badFlies);
    [curvWalks] = removeBadData(curvWalks,badFlies);
    [stopCond] = removeBadData(stopCond,badFlies);
    [boundCond] = removeBadData(boundCond,badFlies);
    
    % again, change the format of the data into empFlys
    [~,~,~,~,~,~,empFlys] = spdCurvProcess(Data,flies,fs,opts.stopThresh,opts.boundary);
    Data = convertSynthetic(empFlys,flies,border);
    [~,speed,~,~,~,~,empFlys] = spdCurvProcess(Data,flies,fs,opts.stopThresh,opts.boundary);
    
    % set thresholds for "run" and "stop" in order to avoid false crossings
    % (note that this is done to be compatable with legacy code from the lab)
    runthreshold = .1;                                                       % 3 mm/s
    stopthreshold = .05;                                                     % 1.5 mm/s
    [~,velocity_classified_binary,~,~] = ...
        Velocity_Classifier({speed'},stopthreshold,runthreshold);
    
    % determine if the fly is inside the light zone or outside the light
    % zone based on head position
    includeZero = false;
    InOuth = {empFlys.rH<border};
    during{1} = false(flies,size(empFlys.rH,2));
    for j = 1:flies
        during{1}(j,empFlys.firstEntry(j):end) = true;
    end
    % get tracks where flies cross the light ring
    [~,Tracks,~,TrackType] = calcCrossings(InOuth,during,velocity_classified_binary,includeZero);
    
    % set the file names to save the data to
    fnameBefore{i} = [newFold '/' filelist{i}(1:end-9) 'CrossingTracksBefore.mat'];
    fnameAfter{i} = [newFold '/' filelist{i}(1:end-9) 'CrossingTracks.mat'];
    % save the crossing data
    saveData(TrackType,Tracks,[1 2],curvPks,curvWalks,stopCond,boundCond,empFlys,fnameBefore{i});
    saveData(TrackType,Tracks,[3 4],curvPks,curvWalks,stopCond,boundCond,empFlys,fnameAfter{i});
    
end

% compute turn bias
saveTurnBias(filelist,fnameBefore,fnameAfter,opts);

end

function [data] = removeBadData(data,badFlies)
% remove any flies that did not enter the center light zone in the during 
% case
data.max(badFlies) = [];
data.tot(badFlies) = [];
data.avg(badFlies) = [];
data.ndx(badFlies) = [];
data.all(badFlies) = [];
data.dirRelCenterBef(badFlies) = [];
data.dirRelCenterAft(badFlies) = [];

end

function [] = saveData(TrackType,Tracks,ndx,curvPks,curvWalks,stopCond,boundCond,empFlys,fileName)
% separate all of the crossing tracks into those that are intering the
% light zone (OutIn) and those that are exiting the light zone (InOUt)

% note that 10798 is hardcoded as that's the number of frames in each track

tmpInOut =TrackType{1,ndx(1)};
tmpTracksInOut = Tracks{1,ndx(1)};
tmpTracksInOut(tmpInOut(:,1)~=ndx(1),:) = [];
tmpInOut(tmpInOut(:,1)~=ndx(1),:) = [];
tmpTracksInOut2 = tmpTracksInOut-tmpInOut(:,3).*10798;

tmpTracksInOut(tmpTracksInOut2>10798) = nan;
tmpTracksInOut2(tmpTracksInOut2>10798) = nan;

tmpTracksInOut = tmpTracksInOut(:,any(~isnan(tmpTracksInOut)));
tmpTracksInOut2 = tmpTracksInOut2(:,any(~isnan(tmpTracksInOut2)));

tmpOutIn =TrackType{1,ndx(2)};
tmpTracksOutIn = Tracks{1,ndx(2)};
tmpTracksOutIn(tmpOutIn(:,1)~=ndx(2),:) = [];
tmpOutIn(tmpOutIn(:,1)~=ndx(2),:) = [];
tmpTracksOutIn2 = tmpTracksOutIn-tmpOutIn(:,3).*10798;

tmpTracksOutIn(tmpTracksOutIn2>10798) = nan;
tmpTracksOutIn2(tmpTracksOutIn2>10798) = nan;

tmpTracksOutIn = tmpTracksOutIn(:,any(~isnan(tmpTracksOutIn)));
tmpTracksOutIn2 = tmpTracksOutIn2(:,any(~isnan(tmpTracksOutIn2)));

save(fileName,'tmpInOut','tmpTracksInOut','tmpTracksInOut2','tmpOutIn',...
    'tmpTracksOutIn','tmpTracksOutIn2','curvPks','curvWalks','stopCond','boundCond','empFlys')
end

function [] = saveTurnBias(filelist,fnameBefore,fnameAfter,opts)
% loop through each crossing file
nGen = length(filelist);
for i = 1:nGen
    % get the turn bias and border choice
    [durProbDetrend,durProb,turnInDuringProb,ratOutIn,ratInOut,turnAmount] = RadDecision(fnameAfter{i},fnameBefore{i},opts);
    save([opts.TBFold '/' filelist{i}(1:end-9) 'BorderChoiceAndTurnBias.mat'],...
        'durProbDetrend','durProb','turnInDuringProb','ratOutIn','ratInOut','turnAmount');
    
    % print figures
    if opts.plotFig
        for k = 1:5
            figure(k)
            suptitle(filelist{i}(1:end-4))
            print('-painters','-dpsc2','AllFigures.ps','-loose','-append');
        end
    end
    close all
end
end






