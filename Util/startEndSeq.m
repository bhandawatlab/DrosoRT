function [startNdx,endNdx,type] = startEndSeq(vec)
startNdx = unique([1 find(diff(vec)~=0)+1]);
endNdx = unique([startNdx(2:end)-1 length(vec)]);
type = vec(startNdx);

% check if the startNdx and endNdx all have the same number of elements
assert(length(startNdx)==length(endNdx));
% check if the endNdx is always after the startNdx
assert(all(endNdx-startNdx>=0));
% check if all of the data is accounted for
assert(sum(endNdx-startNdx+1)==length(vec))
end