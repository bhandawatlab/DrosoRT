function [xErr,yErr] = TSYtoXY(s,fs,sArena)
% This function gets x and y position based on thrust, slip, and yaw
%
% Inputs:
%    s: Structure with fields:
%   `   Kinematics, which in turn has subfields: thrust, slip, and yaw
%   `   Center, which in turn has subfields: x, y (positions)
%    fs: sampling frequency
%    sArena: Structure with field cF to convert from pixels to mm

% Outputs:
%    xErr:
%    yErr: y position discrepency between x,y position and calculated
%    positions

% get thrust,slip, yaw
thrust = s.Kinematics.thrust;
slip = s.Kinematics.slip;
yaw = s.Kinematics.yaw;

% get initial conditions
initAng = s.AngVec(1);
initX = s.Center.x(1);
initY = s.Center.y(1);

% get all orientation in time
allAng = initAng+[0 cumsum(yaw)];
allAng = allAng-360*floor(allAng./360);

% calculate x and y projections of thrust
xThrust = cosd(allAng(1:end-1)).*thrust./fs./sArena.cF;
yThrust = sind(allAng(1:end-1)).*thrust./fs./sArena.cF;

% calculate x and y projections of slip
xSlip = cosd(allAng(1:end-1)-90).*slip./fs./sArena.cF;
ySlip = sind(allAng(1:end-1)-90).*slip./fs./sArena.cF;

% get change in position
delX = xThrust+xSlip;
delY = yThrust+ySlip;

% add delta position to past position
xRec = initX+[0 cumsum(delX)];
yRec = initY+[0 cumsum(delY)];

xAll = s.Center.x;
yAll = s.Center.y;

% see if there are any errors with the empirical x,y
xErr = sum(abs(xAll-xRec));
yErr = sum(abs(yAll-yRec));


end