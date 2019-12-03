function a = myatan(x,y,varargin)
%myatan Choose from various representations of getting the inverse tangent
%
% This function takes in two vectors of dx and dy and returns the vector of
% angles between dx and dy 

% Inputs:
%    Type: 'degrees' or 'radians'
%    Options: 1 = -180:180 or 2 = 0:360

% Test Cases:
% x = [1 -1 0 2 -1 -1 0 -1 sqrt(3)];
% y = [-1 0 1 2 1 -1 -1 -sqrt(3) -1];
% myatan(x,y,'degrees',2);

% 14/08/17 - Initial
% 15/08/17 - added types and options
%
% 2017, Liangyu Tao

if nargin < 2
    error('Must have at least 2 inputs (x,y)')
end
if nargin==2                                                                %just in case the user doesn't give options
    type = 'radians';
    options = 1;
end
nVarargs = length(varargin);
if nVarargs>0 && nVarargs<3
    num = cellfun(@(x) isnumeric(x) && numel(x)==1, varargin);              % true for elements of C that are numerical scalars
    options = [varargin{num},[]];
    type = [varargin{~num},[]];
end

a=zeros(1,length(y));
[~,tmp1] = find(x>0);
[~,tmp2] = find(x==0);
[~,tmp3] = find(x<0);
[~,tmp4] = find(y>=0);
[~,tmp5] = find(y>0);
[~,tmp6] = find(y<0);

case1 = tmp1;
case2 = tmp3(ismember(tmp3,tmp4));
case3 = tmp3(ismember(tmp3,tmp6));
case4 = tmp2(ismember(tmp2,tmp5));
case5 = tmp2(ismember(tmp2,tmp6));

a(case1)=atan(y(case1)./x(case1));
a(case2)=pi+atan(y(case2)./x(case2));
a(case3)=-pi+atan(y(case3)./x(case3));
a(case4)=pi/2;
a(case5)=-pi/2;

if isequal(options,2)
    a(a<0) = a(a<0)+2*pi;
end
if strcmpi(type,'degrees')
    a = a*180/pi;
end

end