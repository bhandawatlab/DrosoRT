function [Crossing,Tracks,TracksBefore,TrackType] = calcCrossings(in,during,velClassBin,includeZero)
% This function finds tracks where the fly is crossing and classifies them
% into crossing in/crossing out and before/during
%
% Inputs:
%    in: cell array where inside each cell is a vector for whether or not
%       the fly is inside the light ring or not
%    during: cell array where inside each cell is a vector for whether or
%       not the light has turned on
%    velClassBin: cell array where inside each cell is a vector for whether
%       or not the fly is moving
%    includeZero: if false, then when the fly is crossing when the fly is
%       defined to be not moving (crossing due to noise), the crossing is
%       discarded
%
% Outputs:
%    Crossing: nx4 cell array where each column of the array contains
%       crossings for a different scenario (before in, before out, during 
%       in, during out) that last timeInt before crossing to 2*timeInt
%       after crossing
%    Tracks: nx4 cell array where each column of the array contains full
%       tracks for a different scenario
%    TracksBefore: nx4 cell array where each column of the array contains
%       a track starting from timeInt before to crossing for the
%       corresponding tracks in Tracks
%    TrackType: for each track, a vector indicating scenario before
%    crossing, scenario after crossing, and when crossing occurs

timeInt = 30;
% calculate crossings of length -timeAfter:timeAfter
for i = 1:length(in)
    numFrames = size(in{i},2);

    crossings = zeros(size(in{i}));
    crossings(in{i} & ~during{i}) = 1;                                      % inside before
    crossings(~in{i} & ~during{i}) = 2;                                     % outside before
    crossings(in{i} & during{i}) = 3;                                       % inside during
    crossings(~in{i} & during{i}) = 4;                                      % outside during
    crossings = crossings';

    crossNdx = find(diff(crossings)~=0);                                    % ndx locations of crossing borders
    
    crossNdx(mod(crossNdx,numFrames)+timeInt>numFrames) = [];               % make sure fly stays in new location for >time after
    crossNdx(mod(crossNdx,numFrames)-timeInt<1) = [];                       % make sure fly was in old location for >time after
    
    
    crossNdx = crossNdx+floor(crossNdx/numFrames);                          % adjust for ndxing
    tempNdx = find(diff(crossNdx)<timeInt);
    crossNdx(tempNdx) = [];                                     % adjust for multiple crossings
    crossNdx(diff(crossNdx)<timeInt) = [];
    
    TimeTrans = repmat(crossNdx,1,timeInt.*2)+[-timeInt+1:1:timeInt];
    
    if includeZero == false
        genBin = velClassBin{i}';genBin = genBin(TimeTrans);
        falseCrossing = genBin(:,timeInt)==0 & genBin(:,timeInt+1)==0;
        TimeTrans(falseCrossing,:) = [];
    end
    
    befScen = crossings(TimeTrans(:,timeInt));
    %aftScen = crossings(TimeTrans(:,timeAfter+1));
    
    Crossing{i,1} = (TimeTrans(befScen==1,:));
    Crossing{i,2} = (TimeTrans(befScen==2,:));
    Crossing{i,3} = (TimeTrans(befScen==3,:));
    Crossing{i,4} = (TimeTrans(befScen==4,:));
    
    In_Out_bef{i} = TimeTrans(befScen==1,:);
    Out_In_bef{i} = TimeTrans(befScen==2,:);
    In_Out_dur{i} = TimeTrans(befScen==3,:);
    Out_In_dur{i} = TimeTrans(befScen==4,:);
end

% calculate all tracks of each scenario
Tracks = cell(length(in),4);TracksBefore = cell(length(in),4);
TrackType = cell(length(in),4);
for i = 1:length(in)
    numFrames = size(in{i},2);

    crossings = zeros(size(in{i}));
    crossings(in{i} & ~during{i}) = 1;                                      % inside before
    crossings(~in{i} & ~during{i}) = 2;                                     % outside before
    crossings(in{i} & during{i}) = 3;                                       % inside during
    crossings(~in{i} & during{i}) = 4;                                      % outside during
    crossings = crossings';
    crossNdx = find(diff(crossings)~=0);                                    % ndx locations of crossing borders
    
    crossNdx(mod(crossNdx,numFrames)+timeInt>numFrames) = [];               % make sure fly stays in new location for >time after
    crossNdx(mod(crossNdx,numFrames)-timeInt<1) = [];                       % make sure fly was in old location for >time after
    crossNdx = crossNdx+floor(crossNdx/numFrames);                          % adjust for ndxing
    tempNdx = find(diff(crossNdx)<timeInt);
    crossNdx(tempNdx) = [];                                     % adjust for multiple crossings
    crossNdx(diff(crossNdx)<timeInt) = [];
    
    
    temp{1} = find(crossings(crossNdx)==1 & crossings(crossNdx+1)==2);
    temp{2} = find(crossings(crossNdx)==2 & crossings(crossNdx+1)==1);
    temp{3} = find(crossings(crossNdx)==3 & crossings(crossNdx+1)==4);
    temp{4} = find(crossings(crossNdx)==4 & crossings(crossNdx+1)==3);
    
    for j = 1:4
        startNdx = crossNdx(temp{j}(temp{j}<length(crossNdx)))+1;
        endNdx = crossNdx(temp{j}(temp{j}<length(crossNdx))+1);
        
        TrackType{i,j} = [crossings(startNdx-1) crossings(startNdx) floor(startNdx/numFrames)];
        Tracks{i,j} = zeros(length(startNdx),max(endNdx-startNdx)+1);
        TracksBefore{i,j} = zeros(length(startNdx),timeInt);
        
        for k = 1:length(startNdx)
            
            currType = crossings(startNdx(k):endNdx(k));
            changeNdx = find(currType ~= currType(1),1);
            if ~isempty(changeNdx)
                endNdx(k) = startNdx(k)+changeNdx-2;
            end
            
            Tracks{i,j}(k,:) = [(startNdx(k):endNdx(k)) NaN(1,size(Tracks{i,j},2)+startNdx(k)-endNdx(k)-1)];
            TracksBefore{i,j}(k,:) = (startNdx(k)-timeInt:startNdx(k)-1);
            
        end
        
    end
end


end