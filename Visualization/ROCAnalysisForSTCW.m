function [] = ROCAnalysisForSTCW(opts)
% This function plots the ROC analysis to determine separability between
% sharp turns and curved walks
%
% Inputs:
%    opts: structure
%       opts.dataConsFold: folder with the consolidated data

close all
% load ST and CW data
if isempty(opts)
    load([pwd '\Data Full\DataCons\Orco Retinal_Nov13.mat'],'curvPks','curvWalks')
else
    load([opts.dataConsFold '\Orco Retinal_Nov13.mat'],'curvPks','curvWalks')
end

predictorTypes = [{'tot Curv'}, {'avg Curv'}];

%obtain total and average curvature (in degrees) for ST and CW
totSTCurv = abs(cell2mat(curvPks.tot)).*180./pi;
totCWCurv = abs(cell2mat(curvWalks.tot)).*180./pi;
avgSTCurv = abs(cell2mat(curvPks.avg)).*180./pi;
avgCWCurv = abs(cell2mat(curvWalks.avg)).*180./pi;

% find the minimum number of data points for the two labels
nSTPt = length(avgSTCurv);
nCWPt = length(avgCWCurv);
nResp = min([nSTPt, nCWPt]);

% sample without replacement such that we have the same number of ST/CW
% labels
avgSTCurv = datasample(avgSTCurv,nResp,'Replace',false);
avgCWCurv = datasample(avgCWCurv,nResp,'Replace',false);
totSTCurv = datasample(totSTCurv,nResp,'Replace',false);
totCWCurv = datasample(totCWCurv,nResp,'Replace',false);

% consolidate the total and avg curv predictors into one label
STPred = [totSTCurv', avgSTCurv'];
CWPred = [totCWCurv', avgCWCurv'];
allPred = [STPred;CWPred];

% define responses, and labels
resp = [false(nResp,1);true(nResp,1)];
labels = cell(nResp.*2,1);
labels(1:nResp) = {'Sharp Turn'};labels(nResp+1:end) = {'Curved Walk'};

figure;
for i = 1:2
    % define predictors
    pred = allPred(:,i);
    
    % conduct logit binary classification followed by a ROC analysis
    mdl = fitglm(pred,resp,'Distribution','binomial','Link','logit');
    scores = mdl.Fitted.Probability;
    [X,Y,~,AUC] = perfcurve(labels,scores,'Curved Walk');
    
    % plot the CDF of the empirical distributions
    subplot(2,2,i.*2-1);hold on
    cdfplot(pred(1:nResp));
    cdfplot(pred(nResp+1:end));
    title(predictorTypes{i})
    xlabel('Degrees');ylabel('CDF')
    
    % plot the ROC curves and the AUC value
    subplot(2,2,i.*2);
    plot(X,Y);
    text(0.25,0.5,['AUC = ' num2str(AUC)])
    xlabel('False positive rate') 
    ylabel('True positive rate')
    title('ROC for Classification by Logistic Regression')
end
set(gcf,'Position',[1.6667, 41.6667, 638.6667, 599.3333])


end



