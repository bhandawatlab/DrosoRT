function [curvPks,curvWalks,stopCond,boundCond] = genSTCWSynth2(synthFlys,type)
% This function generates sharp turn, curved walk, stop, and boundary
% condition structures based on synthetic modeled flies
%
% Inputs
%    empFlys: structure with fields: x,y
%    type: what conditions the fly is in at each time point
% 
% Outputs:
%    curvPks: structure detailing sharp turn condition
%    curvWalks: structure detailing curved walk condition
%    stopCond: structure detailing stop condition
%    boundCond: structure detailing boundary condition

type = type(:,1:end-1);
nFly = size(type,1);
for i = 1:nFly
    s.Center.x = synthFlys.x(i,:);s.Center.y = synthFlys.y(i,:);
    [curv] = CalcCurvature(s);
    [startNdx,endNdx,type2] = startEndSeq(type(i,:));
    
    for j = 1:4
        tmpStart = startNdx(type2==j);tmpEnd = endNdx(type2==j);
        tmpDur = tmpEnd-tmpStart+1;
        
        tmpCurv = zeros(1,length(tmpStart));
        tmpAll = cell(1,length(tmpStart));
        tmpMax = zeros(1,length(tmpStart));
        for k = 1:length(tmpStart)
            a = curv(tmpStart(k):tmpEnd(k));
            tmpCurv(k) = sum(a);
            tmpAll(k) = {[a;tmpStart(k):tmpEnd(k)]'};
            
            [~,ii] = max(abs(a));
            tmpMax(k) = a(ii);
        end
        
        tmpCond{j}.max{i} = tmpMax;
        tmpCond{j}.tot{i} = tmpCurv;
        tmpCond{j}.dur{i} = tmpDur;
        tmpCond{j}.avg{i} = tmpCurv./tmpDur;
        tmpCond{j}.ndx{i} = floor((tmpStart+tmpEnd)./2);
        tmpCond{j}.all{i} = tmpAll;
    end
end

stopCond = tmpCond{1};
curvPks = tmpCond{2};
curvWalks = tmpCond{3};
boundCond = tmpCond{4};

tmp = cell2mat(curvWalks.dur);
tmp2 = cell2mat(curvWalks.ndx);
figure;histogram(tmp(tmp2<1000))

tmp = cell2mat(curvPks.tot).*180./pi;
tmp2 = cell2mat(curvPks.max).*180./pi;
figure;histogram(tmp,[-200:2:200])

tmp = cell2mat(curvWalks.tot).*180./pi;
figure;histogram(tmp(abs(tmp)<200& tmp~=0),[-200:2:200])

tmp = cell2mat(stopCond.tot).*180./pi;
figure;histogram(tmp(abs(tmp)<200& tmp~=0))


end