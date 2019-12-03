function x = discretesample2(p, n)
% process p
K = numel(p);
if ~isequal(size(p), [1, K])
    p = reshape(p, [1, K]);
end

[~,x] = histc(rand(1,n),[0;cumsum(p(:))/sum(p)]);