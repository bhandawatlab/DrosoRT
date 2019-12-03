function [Data] = convertSynthetic(synthFlys,totFly,border)
% This function calculates the curvature and moves the data from a
% structure to a cell array
%
% Inputs:
%    synthFlys: Data structure
%       synthFlys.x, synthFlys.y, synthFlys.r = fly's centroid positions 
%           (n Fly x m time points)
%       synthFlys.xH, synthFlys.yH, synthFlys.rH = fly's head positions  
%           (n Fly x m time points)
%       synthFlys.firstEntry = time index where flies first entered light 
%           ring (n Fly x 1)
%    totFly: total number of flies in synthFlys
%    border: border location

% Outputs:
%    Data: Cell array of data (1 x n Fly)
%       for each cell, in descending rows:
%           x;y;r;in;during;thrust;slip;yaw;curvFly;xH;yH;rH

% calculate the curvature
[curv] = CalcSynCurv(synthFlys,totFly);
Data = cell(1,totFly);
len = size(curv,2);
for fly = 1:totFly
    x = synthFlys.x(fly,1:len);
    y = synthFlys.y(fly,1:len);
    r = synthFlys.r(fly,1:len);
    xH = synthFlys.xH(fly,1:len);
    yH = synthFlys.yH(fly,1:len);
    rH = synthFlys.rH(fly,1:len);
    
    in = rH<border;
    fe = synthFlys.firstEntry(fly);
    during = [false(1,fe-1) true(1,len-fe+1)];
    thrust = sqrt(diff(synthFlys.x(fly,:)).^2+diff(synthFlys.y(fly,:)).^2).*30;
    slip = zeros(1,len);
    yaw = zeros(1,len);
    curvFly = curv(fly,:);
    
    Data{fly} = [x;y;r;in;during;thrust;slip;yaw;curvFly;xH;yH;rH];
end

end

function [curv] = CalcSynCurv(synthFlys,totFly)

curv = zeros(totFly,size(synthFlys.x,2)-1);
for fly = 1:totFly
    x = synthFlys.x(fly,:);
    y = synthFlys.y(fly,:);

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
    curvature = [0 diff(theta)];
    curvature(curvature>pi) = curvature(curvature>pi)-2.*pi;
    curvature(curvature<-pi) = curvature(curvature<-pi)+2.*pi;
    curv(fly,:) = curvature;
end


end
