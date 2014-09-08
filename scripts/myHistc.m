function [h] = myHistc(dat,bins)

% DESCRIPTION:
% Returns an edge-corrected histc
%
% histc finds density  bin(i-1) <= x < bin(i)
% and the last bin counts any values of x that match edges(end). 

% Here, we bundle the number of counts at the edge edges(end) with
% edge(end-1)
% Only the last bin, bin(end-1) <= x <= bin(end)

% HISTORY:
% The need for this function arose in the context of plotting TE and S.
% Strange things were happening at the "zero peak," 
% when AA was zero, making S = 1 (exactly).

% 6/11/12: don't move the density in h(end) to h(end-1)

% REVISION HISTORY:
%
% Ver  Author(s)    Date   Description 
% ---  ----------   ----   --------------------------------------------
%  1    J. Robison  1/4/12  o original code.
%  2    J. Robison  6/8/12  o fixed bug. 
%                           o added comments
%  3    J. Robinson 6/11/12 o don't move the density in h(end) to h(end-1)

% INPUT:
%   dat, the data to be histogramed
%   bins, edges of bins. For 2 bins, have 3 edges....

% OUTPUT:
%   h, the histogrammed data


h = histc(dat,bins);
h(end-1) = h(end-1) + h(end);
h(end) = 0;
return
    
    
    
