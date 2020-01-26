function [] = plotTracksByAttnNdx(opts)
% This function is a wrapper function to plot for different models and
% genotypes, the empirical flies trajectories in order of increasing
% attraction index
%
% Inputs:
%    opts: structure
%       opts.border: light border
%       opts.plotFig: whether or not to plot the figure

%--------------------------------------------------------------------------
% Orco Retinal tracks in ascending fraction time in light zone order
% (Supplementary Figure 6) 
%--------------------------------------------------------------------------
load([pwd '\Run and Tumble\RunMat\Orco Retinal\RT_Kin.mat'],'empFlys')
[attractionNdx] = genAtndx(empFlys,'H',opts);
empFlys = updateFlies(empFlys,attractionNdx);
empFlys.x = empFlys.x./4;
empFlys.y = empFlys.y./4;
fNum = plottingTracks(empFlys,length(empFlys.firstEntry),false,opts,...
    'Orco Retinal (sort by attn ndx)',[]);
if opts.plotFig
    for i = 1:fNum
        figure(i)
        print('-painters','-dpsc2','PaperFigures.ps','-loose','-append');
    end
end

%--------------------------------------------------------------------------
% Orco Control tracks in ascending fraction time in light zone order
%--------------------------------------------------------------------------
load([pwd '\Run and Tumble\RunMat\Orco Control\RT_Kin.mat'],'empFlys')
[attractionNdx] = genAtndx(empFlys,'H',opts);
empFlys = updateFlies(empFlys,attractionNdx);
empFlys.x = empFlys.x./4;
empFlys.y = empFlys.y./4;
fNum = plottingTracks(empFlys,length(empFlys.firstEntry),false,opts,...
    'Orco Control (sort by attn ndx)',[]);
if opts.plotFig
    for i = 1:fNum
        figure(i)
        print('-painters','-dpsc2','PaperFigures.ps','-loose','-append');
    end
end

%--------------------------------------------------------------------------
% Single Antenna tracks in ascending fraction time in light zone order
%--------------------------------------------------------------------------
load([pwd '\Run and Tumble\RunMat\Single Antenna Orco\RT_Kin.mat'],'empFlys')
[attractionNdx] = genAtndx(empFlys,'H',opts);
empFlys = updateFlies(empFlys,attractionNdx);
empFlys.x = empFlys.x./4;
empFlys.y = empFlys.y./4;
fNum = plottingTracks(empFlys,length(empFlys.firstEntry),false,opts,...
    'Single Antennae (sort by attn ndx)',[]);
if opts.plotFig
    for i = 1:fNum
        figure(i)
        print('-painters','-dpsc2','PaperFigures.ps','-loose','-append');
    end
end

%--------------------------------------------------------------------------
% Kinematic model tracks in ascending fraction time in light zone order
% (Supplementary Figure 7) 
%--------------------------------------------------------------------------
load([pwd '\Run and Tumble\RunMat\Orco Retinal\RT_Kin.mat'],'synthFlys')
[attractionNdx] = genAtndx(synthFlys,'H',opts);
synthFlys = updateFlies(synthFlys,attractionNdx);
synthFlys.x = synthFlys.x./4;
synthFlys.y = synthFlys.y./4;
fNum = plottingTracks(synthFlys,length(synthFlys.firstEntry),false,opts,...
    'Kinematic Model (sort by attn ndx)',[]);
if opts.plotFig
    for i = 1:fNum
        figure(i)
        print('-painters','-dpsc2','PaperFigures.ps','-loose','-append');
    end
end

%--------------------------------------------------------------------------
% Kinematic + border choice + turn bias model tracks in ascending fraction 
% time in light zone order (Supplementary Figure 8)
%--------------------------------------------------------------------------
load([pwd '\Run and Tumble\RunMat\Orco Retinal\RT_Kin_BC_TB.mat'],'synthFlys')
[attractionNdx] = genAtndx(synthFlys,'H',opts);
synthFlys = updateFlies(synthFlys,attractionNdx);
synthFlys.x = synthFlys.x./4;
synthFlys.y = synthFlys.y./4;
fNum = plottingTracks(synthFlys,length(synthFlys.firstEntry),false,opts,...
    'Kin+BC+TB Model (sort by attn ndx)',[]);
if opts.plotFig
    for i = 1:fNum
        figure(i)
        print('-painters','-dpsc2','PaperFigures.ps','-loose','-append');
    end
end

end

function [attractionNdx] = genAtndx(DataAll,lab,opts)

r = sqrt(DataAll.(['x' lab]).^2+DataAll.(['y' lab]).^2);
in = r<opts.border*4;
during = false(size(in));
nFlys = size(r,1);
for j = 1:nFlys
    during(j,DataAll.firstEntry(j):end) = true;
end
duringIn = during & in;
beforeIn = ~during & in;
attractionNdx = sum(duringIn,2)./(10800-DataAll.firstEntry);

end

function [synthFlys] = updateFlies(synthFlys,attractionNdx)
[~,ndx] = sort(attractionNdx,'ascend');
synthFlys.x = synthFlys.x(ndx,:);
synthFlys.xH = synthFlys.xH(ndx,:);
synthFlys.y = synthFlys.y(ndx,:);
synthFlys.yH = synthFlys.yH(ndx,:);
synthFlys.r = synthFlys.r(ndx,:);
synthFlys.rH = synthFlys.rH(ndx,:);
synthFlys.thrust = synthFlys.thrust(ndx,:);
synthFlys.slip = synthFlys.slip(ndx,:);
synthFlys.firstEntry = synthFlys.firstEntry(ndx,:);
synthFlys.stopsBefore = synthFlys.stopsBefore(ndx);
synthFlys.stopsDuring = synthFlys.stopsDuring(ndx);
end



