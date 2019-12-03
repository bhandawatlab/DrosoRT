function [data] = normalizeDat(data,n)
% normalizes the data in data to the radius n

data.x = data.x./n;
data.y = data.y./n;
data.r = sqrt(data.x.^2+data.y.^2);

data.xH = data.xH./n;
data.yH = data.yH./n;
data.rH = sqrt(data.xH.^2+data.yH.^2);
end