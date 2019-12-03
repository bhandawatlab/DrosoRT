function [CurvTrack] = computeTurnDirections(x,y,bef,aft,CurvNdx,CurvTrack)
% This function calculates the turning direction
%
% Inputs:
%    x: x-position
%    y: y-position
%    bef: index before the turn
%    aft: index after the turn
%    CurvNdx: index of the turn
%    CurvTrack: structure to save the direction relative to center before
%       and after of turn to
%
% Outputs:
%    CurvTrack: structure to save the direction relative to center before
%       and after of turn to


% compute if the fly is facing more or less towards the center of the arena
u = [x(CurvNdx)-x(bef);y(CurvNdx)-y(bef);zeros(1,length(bef))];
v = [x(aft)-x(CurvNdx);y(aft)-y(CurvNdx);zeros(1,length(bef))];
w = [-x(CurvNdx);-y(CurvNdx);zeros(1,length(bef))];
angBefWalk = zeros(1,size(u,2));angAftWalk = zeros(1,size(u,2));
% calculate the angles between the center facing vector and the
% direction the fly is moving in before the turn and after the turn
for j = 1:size(u,2)
    angBefWalk(j) = atan2(norm(cross(u(:,j),w(:,j))),dot(u(:,j),w(:,j))).*180./pi;
    angAftWalk(j) = atan2(norm(cross(v(:,j),w(:,j))),dot(v(:,j),w(:,j))).*180./pi;
end
CurvTrack.dirRelCenterBef = angBefWalk;
CurvTrack.dirRelCenterAft = angAftWalk;
%--------------------------------------------------------------------------

end