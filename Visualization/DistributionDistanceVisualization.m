function [] = DistributionDistanceVisualization(stopWalkPDFCont,stopWalkPDFRet,turnPDFCont,turnPDFRet,fName)
% This function plots the distributions of each of the 4 kinematic states
% for before, during inside, and during outside
%
% Inputs:
%    stopWalkPDFCont: structure of pdfs for control stop, curved walk, and
%       boundary scenario
%    stopWalkPDFRet: structure of pdfs for retinal stop, curved walk, and
%       boundary scenario
%    turnPDFCont: structure of pdfs for control sharp turn
%    turnPDFRet: structure of pdfs for retinal sharp turn
%    fName: file to print the figures to
%

close all
figNum = 1;subplotNdx = 1;
%d = zeros(8,1);
scenario = {'stop bef','stop out','stop in','before','during out','during in','bound bef','bound dur'};
for i = 1:8
    [nBins,nD] = size(stopWalkPDFCont.val{i});
    dim = nthroot(nBins,nD).*ones(1,nD);
    
    if nD>1
        p1 = reshape(stopWalkPDFCont.prob{i},dim)+eps;
        p2 = reshape(stopWalkPDFRet.prob{i},dim)+eps;
        
        if nD==3
            figNum = figNum+1;subplotNdx = 1;
        else
            if i<4
                figNum = 99;
            else
                figNum = 103;
            end
            subplotNdx = subplotNdx+2;
        end
    else
        figNum = 103;
        p1 = stopWalkPDFCont.prob{i}+eps;
        p2 = stopWalkPDFRet.prob{i}+eps;
        subplotNdx = subplotNdx+1;
    end
    
    p1=p1./sum(p1(:));
    p2=p2./sum(p2(:));
    
    %d(i)=jensen_shannon_divergence(p1,p2);
    
    val = stopWalkPDFCont.val{i};
    val2 = stopWalkPDFRet.val{i};
    plotProbMassFunc(p1,p2,val,val2,nD,dim,figNum,subplotNdx,scenario{i},fName);
end
for i = 99:103
    figure(i)
    if ~isempty(fName)
        print('-painters','-dpsc2',[fName '.ps'],'-loose','-append');
    end
end
% figure(100)
% if ~isempty(fName)
%     print('-painters','-dpsc2',[fName '.ps'],'-loose','-append');
% end

subplotNdx = 3;
for i = 4:6
%     val = turnPDFCont.val{i};
%     val2 = turnPDFRet.val{i};
%     p1 = turnPDFCont.prob{i}./sum(turnPDFCont.prob{i});
%     p2 = turnPDFRet.prob{i}./sum(turnPDFRet.prob{i});
%     
%     subplot(3,2,i);
%     plot(val,p1,'r');hold on;
%     plot(val2,p2,'k')
%     title(scenario{i});
%     xlabel('Curvature')
%     ylabel('Probability')
%     legend('Control','Retinal')
    
    [nBins,nD] = size(turnPDFCont.val{i});
    dim = nthroot(nBins,nD).*ones(1,nD);
    
    p1 = reshape(turnPDFCont.prob{i},dim)+eps;
    p2 = reshape(turnPDFRet.prob{i},dim)+eps;
    p1=p1./sum(p1(:));
    p2=p2./sum(p2(:));
    val = turnPDFCont.val{i};
    val2 = turnPDFRet.val{i};
    
    plotProbMassFunc(p1,p2,val,val2,nD,dim,104,subplotNdx,scenario{i},fName)
%     if i<4
%         xlim([0 360])
%     else
%         xlim([0 180])
%     end
    subplotNdx = subplotNdx+2;
end
suptitle('Sharp Turns')
%set(gcf,'position',[849 49 824 918])
set(gcf,'Position',[1, 41, 638, 599])
if ~isempty(fName)
    print('-painters','-dpsc2',[fName '.ps'],'-loose','-append');
end

end

function [] = plotProbMassFunc(p1,p2,val,val2,nD,dim,figNum,subplotNdx,scenario,fName)

figure(figNum);
yticklabels = [];
if nD == 1
    h = subplot(3,1,subplotNdx);
    plot(val,p1,'r');hold on;
    plot(val2,p2,'k')
    xlim([0 min([find(cumsum(p1)>0.99,1) find(cumsum(p2)>0.99,1)])])
    %set(gcf,'position',[849 49 824 918])
    set(gcf,'Position',[1, 41, 638, 599])
    title(scenario);
    xlabel('Duration')
    ylabel('Probability')
    legend('Control','Retinal')
elseif nD == 2
    subplotNdx = subplotNdx-2;
    h(1) = subplot(3,2,subplotNdx);
    imagesc(p1);
    set(gca, 'YDir', 'Normal');
    title([scenario ' Control'])
    xlabel('Dur (frames)');ylabel('Curv (degrees)');
    
    h(2) = subplot(3,2,subplotNdx+1);
    imagesc(p2);
    set(gca, 'YDir', 'Normal');
    title([scenario ' Retinal'])
    xlabel('Dur (frames)');ylabel('Curv (degrees)');
    
    xticklabels{1} = linspace(min(val(:,1)),max(val(:,1)),dim(1)./10);
    yticklabels{1} = linspace(min(val(:,2)),max(val(:,2)),dim(2)./10);
    xticklabels{2} = linspace(min(val2(:,1)),max(val2(:,1)),dim(1)./10);
    yticklabels{2} = linspace(min(val2(:,2)),max(val2(:,2)),dim(2)./10);
    %     xticklabels{1} = linspace(0,1000,dim(1)./10);
    %     yticklabels{1} = linspace(0,360,dim(2)./10);
    %     xticklabels{2} = linspace(0,1000,dim(1)./10);
    %     yticklabels{2} = linspace(0,360,dim(2)./10);
    %set(gcf,'position',[849 49 824 918])
    set(gcf,'Position',[1, 41, 638, 599])
elseif nD == 3
    subplot(2,1,1)
    plotIso(p1,[scenario ' Control'])
    subplot(2,1,2)
    plotIso(p2,[scenario ' Retinal'])
end

if ~isempty(yticklabels)
    for i = 1:length(h)
        xticks = linspace(1, dim(1), dim(1)./10);
        set(h(i), 'XTick', xticks, 'XTickLabel', xticklabels{i})
        
        yticks = linspace(1, dim(2), dim(2)./10);
        set(h(i), 'YTick', yticks, 'YTickLabel', yticklabels{i})
    end
end

end

function [] = plotIso(pdf,scenario)
ng=75;
MAX = [1,400,10];
MIN = [0.01 5 0];
scaling=MAX-MIN;
% create meshgrid in 3-dimensions
[X1,X2,X3]=meshgrid(MIN(1):scaling(1)/(ng-1):MAX(1),...
    MIN(2):scaling(2)/(ng-1):MAX(2),MIN(3):scaling(3)/(ng-1):MAX(3));

pdf=reshape(pdf,size(X1)); % reshape pdf for use with meshgrid
[B,~] = sort(pdf(:),'descend');
contours_to_consider = {0.5,0.65,0.8,0.9};
alphaVal = 0.7:-0.15:0.25;
cutoff = zeros(size(contours_to_consider));
for j = 1:length(contours_to_consider)
    cutoff(j) = B(find(cumsum(B)>contours_to_consider{j},1));
end

% plot iso surface contours
for j=1:length(cutoff) % isosurfaces with pdf = 0.005,0.01,0.015
    iso = cutoff(j);
    isosurface(X1,X2,X3,pdf,iso),view(3),alpha(alphaVal(j)),box on,hold on
    colormap cool
end
xlabel('Speed (mm/frame)');
ylabel('Duration (frames)');
zlabel('Curvature (deg/frame)');
legend([cellfun(@num2str, contours_to_consider, 'UniformOutput', false) {'data'}]);
title(['Curved walks, ' scenario])
xlim([0 0.75]);ylim([0 100]);zlim([0 3])
view(-150,20)
%set(gcf,'position',[849 49 824 918])
set(gcf,'Position',[1, 41, 638, 599])
end


