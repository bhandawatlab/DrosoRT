addpath(genpath([pwd '\Data Generation']))
addpath(genpath([pwd '\ST Analysis']))
addpath(genpath([pwd '\Visualization']))
addpath(genpath([pwd '\Run and Tumble']))
addpath(genpath([pwd '\Util']))

opts.dataFold = [pwd '\Data Full\Data'];
opts.genFold = [pwd '\Data Full\DataGen'];
opts.STCWFold = [pwd '\Data Full\DataSTCW'];
opts.crossFold = [pwd '\Data Full\DataCrossing'];
opts.TBFold = [pwd '\Data Full\DataTurnBias'];
opts.dataConsFold = [pwd '\Data Full\DataCons'];
opts.border = 1.2/4;% 1.2
opts.stopThresh = 0.5; % in mm/s
opts.boundary = 3.85; % in cm
opts.plotFig = false;

createFolders(opts)

%--------------------------------------------------------------------------
% generate Data Files (new Data) and Data Gen (consolidated Data by genotype)
disp('Generating new data files...')
genGenData(opts);

%--------------------------------------------------------------------------
% % findbest fit parameters to delineate between sharp turn and curved walks
% disp('Finding best parameters...')
% [GlobMinX,GlobMinCFit] = getBestFit(opts);

%--------------------------------------------------------------------------
% generate sharp turn/curved walk data
disp('Generating sharp turn curved walk data files...')
load('BestFit5.mat','GlobMinX');
genSTCWDat(opts,GlobMinX);

%--------------------------------------------------------------------------
% generate crossing and turn bias data
disp('Generating crossing data files...')
genCrossDat(opts);

%--------------------------------------------------------------------------
% if we want to consolidate data
disp('Generating consolidated data files...')
consData(opts);

%--------------------------------------------------------------------------
% Plot basic analysis about distributions for the RT model
file1 = [pwd '\Data Full\DataCons\Orco Control_Nov13.mat'];
file2 = [pwd '\Data Full\DataCons\Orco Retinal_Nov13.mat'];
DistributionExtraction(file1,file2,opts,true);
file1 = [pwd '\Data Full\DataCons\Single Antenna Orco_Nov13.mat'];
file2 = [pwd '\Data Full\DataCons\Orco Retinal_Nov13.mat'];
DistributionExtraction(file1,file2,opts,true);

% Really stupid programming, but I need some functions in DistributionExtraction
% to extract raw values of Orco Kinematics such as speed, curv, etc that is
% separated by state
file = [pwd '\Data Full\DataCons\Orco Retinal_Nov13.mat'];
[~,~,~,rawVals,~,~] = empiricalFlydataGen2(file,opts);
save('Orco Raw Vals.mat', 'rawVals');

%--------------------------------------------------------------------------
% How flies perform ST and comparison to Katsov data
STAnalysisHandle(opts);

--------------------------------------------------------------------------
Create run and tumble model
tic
RunAndTumbleHandle(opts);
toc
%%
%--------------------------------------------------------------------------
% plot figures
%--------------------------------------------------------------------------
opts.plotFig = true;

plotAttractionIndex(opts)
plotChangeInKinematics(opts)
plotDecisionDensity(opts)
plotFirstTurnAnalysis('In2Out',opts)
plotFirstTurnAnalysis('Out2In',opts)
plotSpdCrossing('In2Out',opts)
plotSpdCrossing('Out2In',opts)

files{1} = 'Run and Tumble\RunMat\Orco Retinal\RT_Kin';
files{2} = 'Run and Tumble\RunMat\Orco Retinal\RT_Kin_BC_TB';
opts2.border = opts.border;opts2.gen = 'Orco Retinal';opts2.plotFig = opts.plotFig;
RTPlottingAnalysis(files,opts2);

close all
file1 = [pwd '\Data Full\DataCons\Orco Control_Nov13.mat'];
file2 = [pwd '\Data Full\DataCons\Orco Retinal_Nov13.mat'];
load(file2,'stopWalkPDFs','turnPDFs');
stopWalkPDFs2 =stopWalkPDFs;turnPDFs2 =turnPDFs;
load(file1,'stopWalkPDFs','turnPDFs');
DistributionDistanceVisualization(stopWalkPDFs,stopWalkPDFs2,turnPDFs,turnPDFs2,'PaperFigures')

close all
% plot sample abstracted trajectories (not full tracks)
load('BestFit5.mat', 'GlobMinX')
plotSTCWAbstraction(opts,GlobMinX)
if opts.plotFig
    fNum = get(gcf,'Number');
    for i = 1:fNum
        figure(i);
        print('-painters','-dpsc2','PaperFigures.ps','-loose','-append');
    end
end

% % plot cumsum of curvature
plotCumSumCurv(opts);

% 1 figure for ROC analysis for curved walk and sharp turns
ROCAnalysisForSTCW(opts)
figure(1)
print('-painters','-dpsc2','PaperFigures.ps','-loose','-append');

% 2 figures for a sample track and how to distinguish between sharp turn
% and curved walks
PlotSTCWSegProcess(GlobMinX)
for i = 1:2
    figure(i)
    print('-painters','-dpsc2','PaperFigures.ps','-loose','-append');
end

plotTracksByAttnNdx(opts)

% uncomment this and edit the paths if using ghostscript to convert ps file
% to pdf file
% ps2pdf('psfile', 'PaperFigures.ps', 'pdffile', 'PaperFigures.pdf', 'gspapersize', 'letter',...
% 'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
% 'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
% 'gslibpath','C:\Program Files\gs\gs9.50\lib');


