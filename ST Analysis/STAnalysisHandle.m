function [] = STAnalysisHandle(opts)
% This is the wrapper function for understanding how flies uses thrust,
% slip, and yaw to turn and to compare with the Katsov data
%
% Inputs:
%    opts: structure with fields
%       opts.dataFold: folder with the empirical data by flies
%       opts.STCWFold: folder with sharp turn and curved walk information
%       opts.plotFig: true/false for saving figure plots

yawThresh = 0.0349./2;
curvThresh = 0.0349;
warning('off','signal:findpeaks:largeMinPeakHeight')

% katsov turning
close all
fitFile = 'BestFit5.mat';
[features2,curvPks,curvWalks] = processFeaturesKatsov(yawThresh,curvThresh,fitFile);
spdThresh = [1,10];
SharpTurnAnalysis(features2,spdThresh,[3,2,1]);
spdThresh = [10,20];
SharpTurnAnalysis(features2,spdThresh,[3,2,2]);
spdThresh = [20,100];
SharpTurnAnalysis(features2,spdThresh,[3,2,3]);

% Orco turning
[features] = processFeatures(yawThresh,curvThresh,opts);
% spdThresh = [1,5];
% SharpTurnAnalysis(features,spdThresh,yawThresh,plotFig);
% spdThresh = [5,10];
% SharpTurnAnalysis(features,spdThresh,yawThresh,plotFig);
% spdThresh = [10,100];
% SharpTurnAnalysis(features,spdThresh,yawThresh,plotFig);
spdThresh = [0,1];
SharpTurnAnalysis(features,spdThresh,[3,2,4]);
spdThresh = [1,3];
SharpTurnAnalysis(features,spdThresh,[3,2,5]);
spdThresh = [3,100];
SharpTurnAnalysis(features,spdThresh,[3,2,6]);

compareDistributions(curvPks,curvWalks)

if opts.plotFig
    fNum = get(gcf,'Number');
    for i = 1:fNum
        figure(i);
        print('-painters','-dpsc2','KatsovFigures.ps','-loose','-append');
    end
end

end

