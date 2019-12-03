function [] = checkNChoices(Flys,curvPks,titles)
% This function computes and plots the decision density
%
% Inputs:
%    Flys: cell structure with fields
%       Flys{i}.r: radial location for fly i
%       Flys{i}.firstEntry: time of first entry for fly i
%    curvPks: cell structure with fields
%       curvPks{i}.ndx: time of sharp turns for fly i
%    titles: plotting title
%

% define intensity of light from empirical measurements. For plotting
% purposes
Intensity = [26030, 26020, 26090, 25630, 25800, 25440, 25500, 24920, ...    % Measured intensity in nW
25030, 24780, 25100, 16700, 3610, 134, 31.5, 20];
I = Intensity./max(Intensity);
figure;
for i = 1:length(Flys)
    tmp = [];
    % index where each turn takes place
    for j = 1:length(curvPks{i}.ndx)
        feTmp = Flys{i}.firstEntry(j);
        dTemp = curvPks{i}.ndx{j}(curvPks{i}.ndx{j}-5>=feTmp)-5;
        tmp = [tmp Flys{i}.r(j,dTemp)];
    end
    dx = [0:0.025:1];
    
    % compute radial density of turns and plot them
    decisionDens = histcounts(tmp,dx)./(dx(2:end).^2 - dx(1:end-1).^2);
    normDens = decisionDens./sum(decisionDens);
    subplot(ceil(length(Flys)./2),2,i);plot(dx(1:end-1),normDens);hold on;
    plot([0:0.1/4:1.5/4],I./6);
    plot([1/4 1/4],[0 0.2],'k');
    plot([1.3/4 1.3/4],[0 0.2],'k');
    title(titles{i},'interpreter', 'none')
    
    if i == 1
        legend({'Decision Density','Intensity Profile'})
        xlabel('Normalized DIstance');
        ylabel('Decision Density')
    end
    
end
end