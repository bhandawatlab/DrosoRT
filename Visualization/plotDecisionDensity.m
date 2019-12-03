function [] = plotDecisionDensity(opts)
% This function plots the turn bias and decision density for Orco Retinal
% Flies
%
% Inputs:
%    opts.plotFig: true/false

close all

load([pwd '\Data Full\DataCons\Orco Retinal_Nov13.mat'],'durProb','turnInDuringProb')
x = (0.1:0.1:4.1)/4;
% plot turn bias for the first 2 turns and later turns
figure;subplot(2,1,1);plot(x(1:2:end),turnInDuringProb(1,:),'k');hold on;
plot(x(1:2:end),turnInDuringProb(2,:),'--','Color',[0.5 0.5 0.5]);
plot([1.3/4 1.3/4],[0 1],'r');xlim([0 1]);
ylabel('Probability');xlabel('r Distance');title('Orco Retinal')
legend({'first 2 turns','later turns'})

% plot the probability of turning (normalize for inside) and for outside
% separately
subplot(2,1,2);plot([1:12]./40,durProb(2,1:12)./sum(durProb(2,1:12)),...
    '--','Color',[0.5 0.5 0.5]);hold on;
plot([13:40]./40,durProb(2,13:40)./sum(durProb(2,13:40)),...
    '--','Color',[0.5 0.5 0.5]);
plot([1:12]./40,durProb(1,1:12)./sum(durProb(1,1:12)),'k');hold on;
plot([13:40]./40,durProb(1,13:40)./sum(durProb(1,13:40)),'k');
ylim([0 0.4])
ylabel('Probability');xlabel('r Distance');title('Orco Retinal')
legend({'later turns','later turns','first 2 turns','first 2 turns'})

if opts.plotFig
    n = get(gcf,'Number');
    for i = 1:n
        figure(i);set(gcf,'Position',[2 42 798 774])
        print('-painters','-dpsc2','PaperFigures.ps','-loose','-append');
    end
end

end