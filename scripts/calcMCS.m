function [MCS, photonMap, binMap] = calcMCS(photonArrivalTimes,BIN_TIME,maxTime)

% DESCRIPTION:
% returns the MCS trace from the photon arrival times in global time using
% BIN_TIME binning (usually 1 ms). Also returens the photon assignment maps


% REVISION HISTORY:
%
% Ver  Author(s)    Date   Description 
% ---  ----------   ----   --------------------------------------------
%  1    J. Robison  1/4/12  o formerly getMCS, uses unfiltered events (i.e.
%                            events not yet assigned to  DD, DA, or AA)
%  2    J. Robison  6/8/12  o filtering done separately.
%                           o return map of bin assignments (photonMap and
%                             binMap)


% INPUT:
% photonArrivalTimes, photon arrival times (in global time)
% BIN_TIME, bin time (in msec), usually 1 ms.
% maxTime, ???

% OUTPUT:
%  MCS, the MCS trace
%  photonMap, array of bin assignments
%  binMap, a list of arrays. Each array contains the event assigned to bin.
% binMap = { []
%            [photons x, y, z]      % bin 2
%            []                     % bin 3
%            []
%            ...
%            [] }

numBins = floor(maxTime/BIN_TIME);
BIN_TIME_INC = BIN_TIME * 1E6; % this gives it to us in nsec.
MCS = zeros(numBins,1);
N = length(photonArrivalTimes);
photonMap = zeros(N,1,'uint16');
binMap = cell(numBins,1); 

b = 1; 
counter = 0;
nextTime = BIN_TIME_INC;
myPhotons = [];
for i=1:N
    if photonArrivalTimes(i) < nextTime
        counter = counter + 1; % new photon for bin b                     
        
        % >> append photon i to myPhotons array.
        myPhotons(counter) = i; % bin b has been assigned photon i
    else
        % we've landed in a new interval
        % >> save off the array of positions
        if (length(myPhotons) > 0) % only impossible on the first burst: nothing to save
            binMap{b} = myPhotons;
            myPhotons = []; % register the photon.
            MCS(b) = counter;
        end
        % well, what interval did we land in? We could have an MCS bin
        % with zero counts
        inc = ceil((photonArrivalTimes(i)-nextTime)/BIN_TIME_INC);
        % now get ready for next round
        b = b + inc; % this is the bin we are in.
        counter = 1;% this is the first photon in a new interval
        myPhotons(1) = i;
        nextTime = nextTime + inc*BIN_TIME_INC;
    end;
    photonMap(i) = b; % assign photon i to bin b (which has been updated, if needed)
end;