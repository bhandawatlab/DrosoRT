function [curvature,theta] = CalcCurvatureKatsov(x,y)
% Calculate the curvature of tracks
%
% Inputs:
%    x, y positions

% calculating curvature
Vertices = horzcat(x',y');

%get normal vectors
N=LineNormals2D(Vertices);

theta = zeros(1,length(x)-1); %make it zero vector
for p=1:(length(x)-1)
    %because normal is a unit vector. The x-component can be used
    % to determine the angle it makes with x-axis.
    theta(p)=acos(N(p,1));
    % the if loop converts from 0 to pi to -pi to pi
    if N(p,1)>0 && N(p,2)<0
        theta(p)=-theta(p);
    elseif N(p,1)<0 && N(p,2)<0
        theta(p)=-theta(p);
    end
end

theta(isnan(theta)) = 0; %change Nan to zeros (when the centroid does not move, LineNormals2D assigns 'NaN')
theta(theta<0) = theta(theta<0)+2*pi;
curvature = diff(theta);
curvature = [curvature 0];
curvature(curvature>pi) = curvature(curvature>pi)-2.*pi;
curvature(curvature<-pi) = curvature(curvature<-pi)+2.*pi;


end