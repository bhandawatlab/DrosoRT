function [theta,curv] = norm2Curv(N)

theta = zeros(1,size(N,1)-1); %make it zero vector
for p = 1:length(N)
    theta(p)=acos(N(p,1));
    % the if loop converts from 0 to pi to -pi to pi
    if N(p,1)>0 && N(p,2)<0
        theta(p)=-theta(p);
    elseif N(p,1)<0 && N(p,2)<0
        theta(p)=-theta(p);
    end
end
theta(isnan(theta)) = [];
theta(theta<0) = theta(theta<0)+2*pi;
curv = diff(theta);

curv = [0 curv];
curv(curv>pi) = curv(curv>pi)-2.*pi;
curv(curv<-pi) = curv(curv<-pi)+2.*pi;

end