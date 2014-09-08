% =================================================================
%   TCSPC of PIE-filtered data
% =================================================================

% greenAndRedFilterPos is an array of bins 

green = importdata('../input/APD2_IRF.dat');
red = importdata('../input/IRF of APD1(no filters Sync1).dat');

tcal = 0.016; % 16 ps
sync = 25.0/2; % 25 ns

numPeriods = 6;

chan = floor(sync/tcal);
chan2 = 2*chan;

chans = chan2*numPeriods;
t = tcal*(1:chans);

% assert length(green) < chan

greenTrim = green(1:chan);
redOffset = 770;
redTrim = red((1+redOffset):(chan+redOffset));

figure
hold on
for i = 1:numPeriods
% for i = 1:2
    range = (chan2*(i-1)+1):(i*chan2);
    rangeG = range(1:chan);
    rangeR = range((chan+1):chan2);
    % rangeG([1 end])
    % rangeR([1 end])
    plot(tcal*rangeG,greenTrim,'g','linewidth',1.5);
    plot(tcal*rangeR,redTrim,'r','linewidth',1.5);

end

hold off

set(gca,'FontName','Helvetica');
set(gca,'FontSize',18);
xlabel('Time [nsec]');
ylabel('Intensity (au)');

print -depsc2 '../output/pulseSequence.eps'