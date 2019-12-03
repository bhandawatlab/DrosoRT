function [] = DistributionExtraction(file1,file2,opts,saveDat)
% This function is a wrapper function that generates and saves the 
% probability distributions for each of the 4 states
% 
% Inputs:
%    file1: consolidated data file with the following information
%       'synthFlys','empFlys','curvPks','curvWalks','stopCond','boundCond','gen'
%    file2: consolidated data file with the following information
%       'synthFlys','empFlys','curvPks','curvWalks','stopCond','boundCond','gen'
%    opts: structure with field border for light border
%    saveDat: whether to save distributions to a data file (true/false)
%
[~,stopWalkPDFs,turnPDFs,~,pStop,~] = empiricalFlydataGen2(file1,opts);
[~,stopWalkPDFs2,turnPDFs2,~,pStop2,~] = empiricalFlydataGen2(file2,opts);

% save distributions to the data file
if saveDat
    saveData(file1,stopWalkPDFs,turnPDFs,pStop);
    saveData(file2,stopWalkPDFs2,turnPDFs2,pStop2);
end

end

function [] = saveData(file,stopWalkPDFs,turnPDFs,pStop)
load(file);
save(file);
end


