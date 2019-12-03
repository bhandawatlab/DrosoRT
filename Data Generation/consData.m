function [] = consData(opts)
% This function creates a consolidated file for each genotype with all
% previously calculated information for all flies in that genotype
%
% Inputs:
%    opts: Type of points we want to return (3 options)
%       opts.genFold = folder with the general files of fly movement 
%           information
%       opts.crossFold = folder with the fly crossing light ring
%           information
%       opts.TBFold = folder with the fly turn bias information
%       opts.DataConsFold = folder to save the consolidated files to

% get data folders
DataGenFold = opts.genFold;
DataCrossFold = opts.crossFold;
DataTBFold = opts.TBFold;
DataConsFold = opts.dataConsFold;

% get files in the folders
s = dir([DataGenFold '/*.mat']);
s2 = dir([DataCrossFold '/*.mat']);
s3 = dir([DataTBFold '/*.mat']);

% create cell arrays to enable looping through
filelistGen = {s.name};
filelistCross = {s2.name};
filelistTB = {s3.name};
nGen = length(filelistGen);

% load data from each of the different data folders and assign to a
% consolidated folder
for K = 1:nGen
    load([DataGenFold '/' filelistGen{K}],'Arena','Data','fs','gen');
    load([DataCrossFold '/' filelistCross{K.*2}],'empFlys','curvPks',...
        'curvWalks','stopCond','boundCond');
    load([DataTBFold '/' filelistTB{K}],'turnInDuringProb',...
        'durProbDetrend','durProb','ratInOut','ratOutIn');
    save([DataConsFold '/' filelistGen{K}],'Arena','Data','fs','gen'...
        ,'empFlys','curvPks','curvWalks','stopCond','boundCond',...
        'turnInDuringProb','durProbDetrend','durProb','ratInOut','ratOutIn')
end

end