function MCS = getMCS(photonArrivalTimes,photonFilter, BIN_TIME,maxTime)

% DESCRIPTION:
% returns the MCS trace from the photon arrival times in global time using
% BIN_TIME binning (usually 1 ms).
% photons are selected only if they are within the time gate.

% REVISION HISTORY:
%
% Ver  Author(s)    Date   Description 
% ---  ----------   ----   --------------------------------------------
%  0    P. Kapusta  05/07  o original file name: read_pt3.m (PicoQuant GmbH)
%       M. Wahl
%  1    J. Robison  1/4/12 o Tested, speed up of at least 100x
%  2    J. Robison  6/8/12 o implemented data strcutures.
%                           o writing to a log file.

% INPUT:
% photonArrivalTimes, global photon arrival times
% photonFilter, whether we count a count
% BIN_TIME, in msec

% OUTPUT:
numBins = floor(maxTime/BIN_TIME);
BIN_TIME_INC = BIN_TIME * 1E6; % this gives it to us in nsec.

MCS = zeros(numBins,1);
bin = 1;
counter = 0;
nextTime = BIN_TIME_INC;
for i=1:length(photonArrivalTimes)
    % DD
    if photonFilter(i) % make sure that it passes the time gate
        if photonArrivalTimes(i) < nextTime
            counter = counter + 1;           
        else
            % we've landed in a new interval
            MCS(bin) = counter;
            % well, what interval did we land in? We could have an MCS bin
            % with zero counts
            inc = ceil((photonArrivalTimes(i)-nextTime)/BIN_TIME_INC);
            % now get ready for next round
            bin = bin + inc;
            counter = 1;
            nextTime = nextTime + inc*BIN_TIME_INC;
        end;
    end;
end;