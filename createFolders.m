function [] = createFolders(opts)
% set up what folders to create
fold{1} = opts.genFold;
fold{2} = opts.STCWFold;
fold{3} = opts.crossFold;
fold{4} = opts.TBFold;
fold{5} = opts.dataConsFold;

% create the necessary folder to house the data
for i = 1:numel(fold)
    if ~exist(fold{i}, 'dir')
        mkdir(fold{i})
    end
end

end