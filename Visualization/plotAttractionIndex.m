function [] = plotAttractionIndex(opts)
% This function plots summary statistics for attraction indexes
%
% Inputs:
%    opts: structure with fields
%       opts.border: radial location of light border
%       opts.plotFig: true/false

close all
%--------------------------------------------------------------------------
% Orco Control and Orco Retinal fraction time in light zone (Figure 1c)
%--------------------------------------------------------------------------
file1 = [pwd '\Data Full\DataCons\Orco Control_Nov13.mat'];
file2 = [pwd '\Data Full\DataCons\Orco Retinal_Nov13.mat'];
load(file1,'Data','gen');DataAll{1} = Data;genName{1} = gen.name;
load(file2,'Data','gen');DataAll{2} = Data;genName{2} = gen.name;

[attractionNdx,gLabel] = genLabels(DataAll,'Head',genName,opts);
attractionNdxAll = cell2mat(attractionNdx);
gLabelAll = [gLabel{1};gLabel{2}];
figure;dabest2(attractionNdxAll,gLabelAll,[-0.2 1]);

%--------------------------------------------------------------------------
% Orco Retinal and Single Antenna Orco Retinal fraction time in light zone
%(Figure 6a2)
%--------------------------------------------------------------------------
file1 = [pwd '\Data Full\DataCons\Orco Retinal_Nov13.mat'];
file2 = [pwd '\Data Full\DataCons\Single Antenna Orco_Nov13.mat'];
load(file1,'Data','gen');DataAll{1} = Data;genName{1} = gen.name;
load(file2,'Data','gen');DataAll{2} = Data;genName{2} = gen.name;

[attractionNdx,gLabel] = genLabels(DataAll,'Head',genName,opts);
attractionNdxAll = cell2mat(cellfun(@(x) reshape(x,[],2),attractionNdx,'un',0));
gLabelAll = cellfun(@(x) reshape(x,[],2),gLabel,'un',0);
gLabelAll = [gLabelAll{1};gLabelAll{2}];
gLabelAll = reshape(gLabelAll,[],1);
attractionNdxAll = reshape(attractionNdxAll,[],1);

figure;subplot(2,1,1);
dabest2(attractionNdxAll(1:end/2),gLabelAll(1:end/2),[0 1],'Y');
subplot(2,1,2);
dabest2(attractionNdxAll(end/2+1:end),gLabelAll(end/2+1:end),[-0.2 1],'Y');
set(gcf,'Position',[2 42 798 774])


%--------------------------------------------------------------------------
% Orco Retinal and Synthetic Flies fraction time in light zone
%(Figure 3c, 4c, 5F, Supplement 9)
%--------------------------------------------------------------------------
file{1} = [pwd '\Run and Tumble\RunMat\Orco Retinal\RT_Kin.mat'];
file{2} = [pwd '\Run and Tumble\RunMat\Orco Retinal\RT_Kin_BC_TB.mat'];
for i = 1:3
    if i<3
        load(file{i},'empFlys','synthFlys','gen');
        DataAll{1} = empFlys;genName{1} = 'Orco Retinal Emp';
        DataAll{2} = synthFlys;genName{2} = 'Orco Retinal Synth';
    else
        load(file{1},'synthFlys','gen');
        DataAll{1} = synthFlys;genName{1} = 'Orco Ret Kin Synth';
        load(file{2},'synthFlys','gen');
        DataAll{2} = synthFlys;genName{2} = 'Orco Ret Full Synth';
    end
    DataAll{1}.lightOn = DataAll{1}.firstEntry';
    DataAll{2}.lightOn = DataAll{2}.firstEntry';
    
    [attractionNdx,gLabel] = genLabels(DataAll,'H',genName,opts);
    attractionNdxAll = cell2mat(cellfun(@(x) reshape(x,[],2),attractionNdx,'un',0));
    gLabelAll = cellfun(@(x) reshape(x,[],2),gLabel,'un',0);
    gLabelAll = [gLabelAll{1};gLabelAll{2}];
    gLabelAll = reshape(gLabelAll,[],1);
    attractionNdxAll = reshape(attractionNdxAll,[],1);
    
    if i<3
        figure;subplot(2,1,1);
        dabest2(attractionNdxAll(1:end/2),gLabelAll(1:end/2),[0 1],'Y');
        subplot(2,1,2);
        dabest2(attractionNdxAll(end/2+1:end),gLabelAll(end/2+1:end),[-0.2 1],'Y');
    else
        sce2cons = {'Orco Ret Kin Synth Before','Orco Ret Kin Synth During','Orco Ret Full Synth During'};
        Index= [];
        for i = 1:length(sce2cons)
            Index = [Index; find(contains(gLabelAll,sce2cons{i}))];
        end
        g2Cons = gLabelAll(Index);attractionNdx2Cons = attractionNdxAll(Index);
        figure;
        dabest2(attractionNdx2Cons,g2Cons,[-0.2 1],'N');
    end
    %set(gcf,'Position',[2 42 838 924])
end


% file{1} = [pwd '\Run and Tumble\RunMat\Orco Retinal\RT_Kin.mat'];
% file{2} = [pwd '\Run and Tumble\RunMat\Orco Retinal\RT_Kin_BC.mat'];
% file{3} = [pwd '\Run and Tumble\RunMat\Orco Retinal\RT_Kin_BC_TB.mat'];
% file{4} = [pwd '\Run and Tumble\RunMat\Orco Retinal\RT_Kin_BC_TB.mat'];
% for i = 1:4
%     load(file{i},'synthFlys','empFlys');DataAll{i} = synthFlys;genName{i} = file{i}(85:end-4);
%     DataAll{i}.lightOn = DataAll{i}.firstEntry';
% end
% DataAll{5} = empFlys;genName{5} = 'Orco Retinal';
% DataAll{5}.lightOn = DataAll{5}.firstEntry';
%
% [attractionNdx,gLabel] = genLabels(DataAll,'H',genName,opts);
% attractionNdxAll = cell2mat(cellfun(@(x) reshape(x,[],2),attractionNdx,'un',0));
% gLabelAll = cellfun(@(x) reshape(x,[],2),gLabel,'un',0);
% gLabelAll = [gLabelAll{1};gLabelAll{2};gLabelAll{3};gLabelAll{4};gLabelAll{5}];
% gLabelAll = reshape(gLabelAll,[],1);
% attractionNdxAll = reshape(attractionNdxAll,[],1);
%
% % figure;subplot(2,1,1);
% % dabest2(attractionNdxAll(1:end/2+35),gLabelAll(1:end/2+35),[0 1],'Y');
% figure;
% dabest2(attractionNdxAll(end/2+1-35:end),gLabelAll(end/2+1-35:end),[0 1],'Y');
% set(gcf,'Position',[2 42 838 924])
% print('-painters','-dpsc2','DiscoverDay.ps','-loose','-append');

if opts.plotFig
    n = get(gcf,'Number');
    for i = 1:n
        figure(i);set(gcf,'Position',[2 42 798 774])
        suptitle('Attraction Index (fraction time in light zone)')
        %print('-painters','-dpsc2','PaperFigures.ps','-loose','-append');
    end
end

end



function [attractionNdx,gLabel] = genLabels(DataAll,lab,genName,opts)

attractionNdx = cell(2,1);gLabel = cell(2,1);
for i = 1:numel(DataAll)
    r = sqrt(DataAll{i}.(['x' lab]).^2+DataAll{i}.(['y' lab]).^2);
    in = r<opts.border*4;
    during = false(size(in));
    nFlys = size(r,1);
    for j = 1:nFlys
        during(j,DataAll{i}.lightOn(j):end) = true;
    end
    duringIn = during & in;
    beforeIn = ~during & in;
    attractionNdx{i} = [sum(beforeIn,2)./(10800-DataAll{i}.lightOn');
        sum(duringIn,2)./(10800-DataAll{i}.lightOn')];
    gLabel{i} = [repmat({[genName{i} ' Before']},nFlys,1);
        repmat({[genName{i} ' During']},nFlys,1)];
end

end











