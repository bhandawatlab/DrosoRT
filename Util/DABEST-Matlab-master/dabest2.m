function [ss,avr,moes] = dabest2(data,identifiers,varargin)

if size(data,2) ~=1
    identifiers = repmat(identifiers,size(data,1),1);
    data = reshape(data,[],1);
    identifiers = reshape(identifiers,[],1);
end

avr = [];moes = [];
% if ~isempty(varargin)
%     if strcmp(varargin{1},'Paired')
%         [ss] = FscatJit2(identifiers, data,'Y');
%     
%     elseif strcmp(varargin{1},'mergeGroups')
%         [ss,avr,moes] = FscatJit2_mergeGroups(identifiers, data);
%     end
% else
    [ss] = FscatJit2(identifiers, data,varargin{1});
% end

end