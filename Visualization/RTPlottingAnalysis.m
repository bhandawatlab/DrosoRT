function [flies] = RTPlottingAnalysis(files,opts)
% This function is a wrapper function to generate some simple analysis for
% the run and tumble models and plot radial occupancy and where they make
% turns
%
% Inputs:
%    files: cell array of file names for RT model outputs
%    opts: structure
%       opts.border: light border
%       opts.plotFig: whether or not to plot the figure

close all
nFiles = length(files);
border = opts.border;
%opts.gen = 'Orco Retinal';
% loop through each synthetic model to compute some basic analysis and
% sharp turn/curved walks
for i = 1:nFiles
    [flies{i},~,~,~,~,~] = empiricalFlydataGen2(files{i},opts,'Synth');
    fileName{i} = files{i}(39:end);
    nFlys{i} = length(flies{i}.firstEntry);
    [flyAnalysis{i}] = basicAnalysis(flies{i},nFlys{i},border);
    load(files{i},'curvPks');
    turnPDFs{i} = curvPks;
end
load([pwd '/Data Full/DataCons/' opts.gen '_Nov13'],'curvPks');
turnPDFs{nFiles+1} = curvPks;

% generate sharp turn and curved walk distributions and some basic analysis
% for empirical flies
[flies{nFiles+1},~,~,~,~,~] = empiricalFlydataGen2(files{nFiles},opts,'Emp');
fileName{nFiles+1} = 'Empirical';
nFlys{nFiles+1} = length(flies{nFiles+1}.firstEntry);
[flyAnalysis{nFiles+1}] = basicAnalysis(flies{nFiles+1},nFlys{nFiles+1},border);

% plotting functions
plotRadialOccupancy(flyAnalysis,opts,fileName);
checkNChoices(flies,turnPDFs,fileName);

if opts.plotFig
    n = get(gcf,'Number');
    for i = 1:n
        figure(i);set(gcf,'Position',[2 42 798 774])
        print('-painters','-dpsc2','PaperFigures.ps','-loose','-append');
    end
end

end

function [flies,stopWalkPDFs,turnPDFs,rawVals,pStop,gen] = empiricalFlydataGen2(file,opts,type)
if strcmpi(type,'Emp')
    load(file,'empFlys','curvPks','curvWalks','stopCond','boundCond','gen')
    flies = empFlys;
else
    load(file,'synthFlys','curvPks','curvWalks','stopCond','boundCond','gen')
    flies = synthFlys;
end

if ~exist('gen','var')
    gen = [];
end

% get raw values
[flies,spdDist,durDist,durDist2,yawDist,curvDist,pStop] = ...
    rawValExtraction(flies,curvPks,curvWalks,stopCond,boundCond,opts);

% get probability distributions
stopWalkPDFs = DistDurCurvRelationship2(spdDist,durDist,curvDist);
turnPDFs = turnPDFGen(yawDist,durDist2);

rawVals.spd = spdDist;
rawVals.dur = durDist;
rawVals.sharpTurn = yawDist;
rawVals.curvWalk = curvDist;

end

