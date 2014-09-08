function y = filterSelect(data,sel)
% DESCRIPTION
%   fill _f with members of _data that have _sel = 1. 
% INPUT
%   data, data
%   sel, selection filter
% OUTPUT
%   filter

if length(data) ~= length(sel)
    error('DataVectorSizesDoNotMatch','The number of elements in dta and sel do not match')
end
%f = double(1);

% only export those within the filter
j = 1;
for i = 1:length(data)
    if sel(i)
        f(j)= data(i);
        j = j +1;        
    end
end
y = transpose(f);

