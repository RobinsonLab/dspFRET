function fP = filterPos(filter)

% DESCRIPTION:
% return maps that will be used for photon selection
% output the (sequential) locations of True in an input array


% REVISION HISTORY:
%
% Ver  Author(s)    Date   Description 
% ---  ----------   ----   --------------------------------------------
%  1    J. Robison  6/8/12 o first implementation

% INPUT:
%  filter, boolean array of picked points. 

% OUTPUT:
%  positions of T.

n = length(filter);
c = 0;
m = sum(filter);
fP = zeros(m,1,'uint16'); 
for i = 1:n % strange, but i is double.
    if filter(i)
        c = c + 1;
        fP(c) = i; % will matlab cast i into uint16? 
    end
end
