function [Data,badFlies] = convertFormat(dat,Arena,border,nFlies)
% This function moves fly data from a structure to a cell array
%
% Inputs:
%    dat: Data structure
%       dat.x, dat.y = fly's centroid positions (n Fly x m time points)
%       dat.xHead, dat.yHead = fly's head positions (n Fly x m time points)
%       dat.lightOn = time index when light is turn on
%       dat.thrust, dat.slip, dat.yaw, dat.curv = fly's kinematics 
%           (n Fly x m time points)
%    Arena: struction
%       Arena.rad: radius of arena
%       Arena.cF: conversion factor from pixels to mm
%    border: border location
%    nFlies: number of flies
%
% Outputs:
%    Data: Cell array of data (1 x n Fly)
%       for each cell, in descending rows:
%           x;y;r;in;during;thrust;slip;yaw;curv;xH;yH;rH
%    badflies: indexes of flies that never entered the light ring during
%       the light on period

Data = cell(1,nFlies);
for i = 1:nFlies
    x = dat.x(i,:);
    y = dat.y(i,:);
    xH = dat.xHead(i,:);
    yH = dat.yHead(i,:);
    r = sqrt(x.^2+y.^2);
    rH = sqrt(xH.^2+yH.^2);
    
    in = rH<=border;
    
    fe = find(in(dat.lightOn(i):end),1)+dat.lightOn(i)-1;
    
    during = [false(1,fe-1) true(1,length(in)-fe+1)];
    
    thrust = dat.thrust(i,:);%./30./Arena.rad(i)./Arena.cF(i);
    slip = dat.slip(i,:);%./30./Arena.rad(i)./Arena.cF(i);
    yaw = dat.yaw(i,:);
    curv = dat.curv(i,:);
    
    if ~isempty(during)
        Data{i} = [x(1:end-1);y(1:end-1);r(1:end-1);in(1:end-1);during(1:end-1);...
            thrust;slip;yaw;curv;xH(1:end-1);yH(1:end-1);rH(1:end-1)];
    end
end
badFlies = find(cellfun(@isempty,Data));
Data(cellfun(@isempty,Data)) = [];

end