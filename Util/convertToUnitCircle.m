function [Data] = convertToUnitCircle(Data,Arena)

for i = 1:size(Data.Center.x,1)
    xTmp = Data.Center.x(i,:);
    yTmp = Data.Center.y(i,:);
    x = interp1q(find(xTmp~=0)',xTmp(xTmp~=0)',[1:1:length(xTmp)]')';
    y = interp1q(find(yTmp~=0)',yTmp(yTmp~=0)',[1:1:length(yTmp)]')';
    Data.Center.x(i,:) = x;Data.Center.y(i,:) = y;
end

x = Data.Center.x;y = Data.Center.y;

xH = Data.Head.x;
yH = Data.Head.y;

for i = 1:size(x,1)
    x(i,:) = (x(i,:)-Arena.arenaCent(1,i))./Arena.rad(i);
    y(i,:) = (y(i,:)-Arena.arenaCent(2,i))./Arena.rad(i);
    xH(i,:) = (xH(i,:)-Arena.arenaCent(1,i))./Arena.rad(i);
    yH(i,:) = (yH(i,:)-Arena.arenaCent(2,i))./Arena.rad(i);
    
    rad = sqrt(x(i,:).^2+y(i,:).^2);
    
    if sum(rad>0.85)>100
        x(i,:) = x(i,:)./max(rad);
        y(i,:) = y(i,:)./max(rad);
        xH(i,:) = xH(i,:)./max(rad);
        yH(i,:) = yH(i,:)./max(rad);
    end
    
end

Data.Center.x = x;
Data.Center.y = y;
Data.Head.x = xH;
Data.Head.y = yH;

end
