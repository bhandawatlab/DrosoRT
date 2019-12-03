function [] = RunAndTumbleHandle(fileOpts)
% This function is the handle code to set up and run the different
% complexities of "run and tumble" model
%
% Inputs:
%    fileOpts: Type of points we want to return (3 options)
%       opts.border = normalized radial position of the light border (0>=b<=1)
%       opts.dataConsFold = folder with the consolidated data files

% set border
options.border = fileOpts.border;
% set folder to save the models into
options.folder = 'Run and Tumble';

% "kinematic" Run and Tumble (RT)
options.kin = true;
options.stops = true;
options.curvedWalks = true;
options.borderChoice = false;
options.SharpTurnBiase = false;
options.CurvedWalkBiase = false;
options.noise = false;
opts{1} = options;

% RT with Border Choice
options.borderChoice = true;
opts{2} = options;

% RT with turn bias
options.SharpTurnBiase = true;
options.CurvedWalkBiase = true;
opts{3} = options;

% RT with only turn biases, no border choice
options.borderChoice = false;
opts{4} = options;

% RT with decision space, but not kinematic space
options.borderChoice = true;
options.kin = false;
opts{5} = options;

warning('off');
dataFold =dir(fullfile(fileOpts.dataConsFold, '*.mat'));
dataFiles = {dataFold.name};
% cell array of file name identifiers to save the models to
fName = {'Kin','Kin_BC','Kin_BC_TB','Kin_TB','BC_TB'};
%--------------------------------------------------------------------------
% start a parallel pool. Use this out if not using parallel programing
% delete(gcp('nocreate'))
% parpool('local',6)
%--------------------------------------------------------------------------
% %ppm = ParforProgMon('Processing: ', length(files) , 1);

% loop through each genotype
for i = 1:length(dataFiles)
    inputFile = [fileOpts.dataConsFold '\' dataFiles{i}];
    currGen = dataFiles{i}(1:end-10);
    mkdir([options.folder '\RunMat\' currGen])
    %mkdir([options.folder '\RunMat\' num2str(trial) '\' currGen])
    
    % loop through each RT scenario
    for j = 1:length(opts)
        close all
        options = opts{j};
        options.fileName = [options.folder '\RunMat\' currGen '\RT_' fName{j}];
        %options.fileName = [options.folder '\RunMat\' num2str(trial) '\' currGen '\RT_' fName{j}];
        
        % run the run and tumble model wrapper
        RunAndTumbleWrapper(inputFile,options,fileOpts);
    end
end
toc

% loop through each genotype to generate analysis figures for each model
% fit
for i = 1:length(dataFiles)
    opts2.border = fileOpts.border;
    opts2.gen = dataFiles{i}(1:end-10);
    opts2.plotFig = false;
    files = cell(1,length(fName));
    % get file names
    for j = 1:length(fName)
        files{j} = [options.folder '/RunMat/' opts2.gen '/RT_' fName{j}];
        %files{j} = [options.folder '\RunMat\' num2str(trial) '\' opts2.gen '\RT_' fName{j}];
    end
    % order of increasing complexity (naming convention makes it out of order)
    files = files([1,2,4,5,3]);
    % plot the analysis figures
    RTPlottingAnalysis(files,opts2);
    
    n = get(gcf,'Number');
    for j = 1:n
        figure(j);set(gcf,'Position',[2 42 798 774])
        print('-painters','-dpsc2',[options.folder '\RunMat\' opts2.gen '/PaperFigures.ps'],'-loose','-append');
        %print('-painters','-dpsc2',[options.folder '\RunMat\' num2str(trial) '\' opts2.gen '/PaperFigures.ps'],'-loose','-append');
    end
    
    close all
end
ps2pdf('psfile', [options.folder '\RunMat\' opts2.gen '/PaperFigures.ps'],...
'pdffile', [options.folder '\RunMat\PaperFigures.pdf'], 'gspapersize', 'letter',...
'gscommand','C:\Program Files\gs\gs9.50\bin\gswin64.exe',...
'gsfontpath','C:\Program Files\gs\gs9.50\lib',...
'gslibpath','C:\Program Files\gs\gs9.50\lib');

end





