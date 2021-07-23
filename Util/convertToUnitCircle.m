function [Data] = convertToUnitCircle(Data,Arena)

for i = 1:size(Data.x,1)
    xTmp = Data.x(i,:);
    yTmp = Data.y(i,:);
    x = interp1q(find(xTmp~=0)',xTmp(xTmp~=0)',[1:1:length(xTmp)]')';
    y = interp1q(find(yTmp~=0)',yTmp(yTmp~=0)',[1:1:length(yTmp)]')';
    Data.x(i,:) = x;Data.y(i,:) = y;
end

x = Data.x;y = Data.y;

xH = Data.xHead;
yH = Data.yHead;

for i = 1:size(x,1)
    x(i,:) = (x(i,:)-Arena.center(1,i))./Arena.rad(i);
    y(i,:) = (y(i,:)-Arena.center(2,i))./Arena.rad(i);
    xH(i,:) = (xH(i,:)-Arena.center(1,i))./Arena.rad(i);
    yH(i,:) = (yH(i,:)-Arena.center(2,i))./Arena.rad(i);
    
    rad = sqrt(x(i,:).^2+y(i,:).^2);
    
    if sum(rad>0.9)>100
        x(i,:) = x(i,:)./max(rad);
        y(i,:) = y(i,:)./max(rad);
        xH(i,:) = xH(i,:)./max(rad);
        yH(i,:) = yH(i,:)./max(rad);
    end
    
end

Data.x = x;
Data.y = y;
Data.xHead = xH;
Data.yHead = yH;

end