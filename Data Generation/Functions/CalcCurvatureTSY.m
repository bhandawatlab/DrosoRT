function [allCurv,err,yaw] = CalcCurvatureTSY(s,fs,sArena)
% This function calculates the curvature and the normal vector based on the
% thrust (T), slip (S), and yaw (Y)
%
% Inputs:
%    s: cell array with 
%       position fields: s.Center.x, s.Center.y
%       kinematic fields: s.Kinematics.thrust,
%           s.Kinematics.slip,s.Kinematics.yaw
%    fs: sampling frequency
%    sArena: cell array with field
%       cF: convert from pixels to mm
% 
% Outputs:
%    allCurv: curvature of the trajectory based on x,y positions and TSY
%       separately
%    err: error between the normals calculated by the two methods
%    yaw: s.Kinematics.yaw

x = s.Center.x;
y = s.Center.y;

% calculating curvature
Vertices = horzcat(x',y');

%------------------------------------------------------------------------
% calculate normals using one of two methods (both are the same)
% Method 1
%------------------------------------------------------------------------
N=LineNormals2D(Vertices);
% Method 2
%------------------------------------------------------------------------
dVert = [[0,0]; -diff(Vertices); [0,0]];
N1 = dVert(1:end-1,:)./(sum(dVert(1:end-1,:)'.^2)')+...
    dVert(2:end,:)./(sum(dVert(2:end,:)'.^2)');
N1 = [-N1(:,2) N1(:,1)]./sqrt(sum(N1'.^2)');
dxy = diff(Vertices);
%------------------------------------------------------------------------

% calculate curvature and theta from the normals
[theta,curv] = norm2Curv(N1);

% define thrust, slip, yaw, orientation
thrust = s.Kinematics.thrust./fs./sArena.cF;
slip = s.Kinematics.slip./fs./sArena.cF;
yaw = s.Kinematics.yaw;
initAng = s.AngVec(1);

cumYaw = initAng+[0 cumsum(yaw(1:end-1))];

% calculate dx and dy
dx = thrust.*cosd(cumYaw)+slip.*sind(cumYaw);
dy = thrust.*sind(cumYaw)-slip.*cosd(cumYaw);
s = (dx.^2+dy.^2);

dxs = dx(1:end-1)./s(1:end-1)+dx(2:end)./s(2:end);
dys = dy(1:end-1)./s(1:end-1)+dy(2:end)./s(2:end);

% calculate normals
N2 = [dys',-dxs'];
N2 = N2./sqrt(sum(N2'.^2)');

[theta2,curv2] = norm2Curv(N2);

% double check thetas
mLen = min(length(curv),length(curv2));
err = sum(abs(curv(1:mLen)-curv2(1:mLen)));

allCurv = [curv(1:mLen)',curv2(1:mLen)'];

end