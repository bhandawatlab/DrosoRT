function [] = plotRadialOccupancy(flyAnalysis,opts,fileName)
% This function plots the radial occupancy
%
% Inputs:
%    flyAnalysis: cell structure with fields
%       Flys{i}.Prob.before: x = radial position, y = probability
%       Flys{i}.Prob.during: x = radial position, y = probability
%    opts: structure with fields
%       opts.border: radial location of border
%       opts.gen: genotype name
%    fileName: name of file

border = opts.border;
figure;subplot(2,1,1)
set(gcf,'Position',[1.6667, 41.6667, 638.6667, 599.3333])
c = varycolor(length(flyAnalysis));
for i = 1:length(flyAnalysis)
    % plot radial probability
    plot(flyAnalysis{i}.Prob.before.x,flyAnalysis{i}.Prob.before.y,'Color',c(i,:),'Linewidth',1);hold on;
end
plot([border,border],[0 0.6],'r--');
legend([fileName {'border'}],'Location','northEast','Interpreter', 'none')
xlabel('Normalized radial distance');ylabel('Probability');title('Before Radial Occupancy')

subplot(2,1,2)
for i = 1:length(flyAnalysis)
    % plot radial probability
    plot(flyAnalysis{i}.Prob.during.x,flyAnalysis{i}.Prob.during.y,'Color',c(i,:),'Linewidth',1);hold on;
end
plot([border,border],[0 0.6],'r--');
legend([fileName {'border'}],'Location','northEast','Interpreter', 'none')
xlabel('Normalized radial distance');ylabel('Probability');title('During Radial Occupancy')
suptitle(opts.gen)

end