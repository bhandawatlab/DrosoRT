function [pdfs] = turnPDFGen(yawDist,durDist)
% This function calculates joint probability mass functions for sharp turns
% 
% Inputs:
%    yawDist: 1x6 cell array curvature of stops and sharp turns
%    durDist: 1x6 cell array duration of stops and sharp turns
% * note that in order, these scenarios are:
%       before stop; during Outside stop; during Inside stop; before sharp 
%       turn; during Outside sharp turn; during Inside sharp turn;

% Outputs:
%    pdfs: Structure of the pdfs for each scenario
%       pdfs.type: 1xn cell array to label what each dimension means
%       pdfs.prob: linearized probability density (m grids x n dimensions)
%       pdfs.val: corresponding values for each grid to pdfs.prob

pdfs = [];
% for i = 1:6
%     data = yawDist{i}';
%     if isempty(data)
%         data = 0;
%     end
%     %grid = [min(data):1:max(data)]';
%     grid = [0:1:360]';
%     pd = fitdist(data,'Kernel');
%     pdf2 = pdf(pd,grid);
%     pdfs.val{i} = grid;
%     pdfs.type{i} = {'Turn'};
%     pdf2 = pdf2./sum(pdf2);
%     pdfs.prob{i} = pdf2;
% end

for i = 1:6
    ng = 75;
    if i < 4
        MAX=[1800 360]; MIN=[0 0];
    else
        MAX=[90 180]; MIN=[0 0];
    end
    
    data = [durDist{i}',yawDist{i}'];
    if isempty(data)
        data = [0,0];
    end
    [~,d]=size(data);
    %MAX=max(data,[],1); MIN=min(data,[],1);
    scaling=MAX-MIN;
    [X1,X2]=meshgrid(MIN(1):scaling(1)/(ng-1):MAX(1),...
        MIN(2):scaling(2)/(ng-1):MAX(2));
    grid=reshape([X1(:),X2(:)],ng^d,d);
    pdf2=akde(data,grid); % run adaptive kde
    pdfs.type{i} = {'Dur','Curv'};
    
    pdf2 = pdf2./sum(pdf2);
    pdfs.prob{i} = pdf2;
    pdfs.val{i} = grid;
end


end