function [mm] = meanmin(vec,nMINS)
%MULTIMIN return mean of the lowest nMINS valus in vec

for c=1:nMINS
    [tmp(c),i] = min(vec);
    vec(i) = max(vec);
end
mm=mean(tmp);
