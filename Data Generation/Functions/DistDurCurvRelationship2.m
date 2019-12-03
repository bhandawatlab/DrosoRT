function [pdfs] = DistDurCurvRelationship2(spdDist,durDist,curvDist)
% This function calculates joint probability mass functions for stops,
% curved walks, and boundary conditions
%
% Inputs:
%    spdDist: 1x8 cell array of avg speed for each instance of each
%    scenario
%    durDist: 1x8 cell array of duration for each instance of each
%    scenario
%    curvDist: 1x8 cell array of avg curvature for each instance of each
%    scenario
% * note that in order, these scenarios are:
%       before stop; during Outside stop; during Inside stop; before run; 
%       during Outside run; during Inside run; before boundary; during boundary;

% Outputs:
%    pdfs: Structure of the pdfs for each scenario
%       pdfs.type: 1xn cell array to label what each dimension means
%       pdfs.prob: linearized probability density (m grids x n dimensions)
%       pdfs.val: corresponding values for each grid to pdfs.prob

pdfs = [];
for i = 1:8
    ng=75; % total grid points = ng^d
%     if i<4
%         data = durDist{i}';
%         if isempty(data)
%             data = 0;
%         end
%         grid = [min(data):1:max(data)]';
%         pd = fitdist(data,'Kernel');
%         pdf2 = pdf(pd,grid);
%         pdfs.type{i} = {'Dur'};
%         
    % for curved walks
    if i<7 && i>3
            data = [spdDist{i}',durDist{i}',curvDist{i}'];
            [~,d]=size(data);MAX = [1,700,15];MIN = [0.01 5 0];
            scaling=MAX-MIN;
            % create meshgrid in 3-dimensions
            [X1,X2,X3]=meshgrid(MIN(1):scaling(1)/(ng-1):MAX(1),...
                MIN(2):scaling(2)/(ng-1):MAX(2),MIN(3):scaling(3)/(ng-1):MAX(3));
            grid=reshape([X1(:),X2(:),X3(:)],ng^d,d);
            pdf2=akde(data,grid); % run adaptive kernel density estimation
            pdfs.type{i} = {'Dur','Spd','Curv'};
            grid = grid(:,[2,1,3]);
    else
        % for sharp turns and stop conditions
        data = [durDist{i}',curvDist{i}'];
        if ~isempty(data)
            [~,d]=size(data);
            %MAX=max(data,[],1); MIN=min(data,[],1); 
            MAX=[1800 360]; MIN=[0 0]; 
            scaling=MAX-MIN;
            % create meshgrid in 2 dimensions
            [X1,X2]=meshgrid(MIN(1):scaling(1)/(ng-1):MAX(1),...
                MIN(2):scaling(2)/(ng-1):MAX(2));
            grid=reshape([X1(:),X2(:)],ng^d,d);
            pdf2=akde(data,grid); % run adaptive kde
            pdfs.type{i} = {'Dur','Curv'};
        else
            pdf2 = ones(ng.^2,1)./(ng.^2);
        end
    end
    
    pdf2 = pdf2./sum(pdf2);
    pdfs.prob{i} = pdf2;
    pdfs.val{i} = grid;
end

end