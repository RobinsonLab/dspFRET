function plotBurst1(tstep,dat,start,len,str)
% plots both donor and acceptor channels in a burst trace. 
% those bursts that are selected by the PIE filter are shown with a grey
% background.


% INPUTS:
% tstep, binning time step [msec]
% don, a array of size (t,2) 
% start, start channel to plot data
% len, number of bins to plot
% str, ylabel

% OUTPUTS:



% figure;

% CONSTANTS
fac = 1.2;
c = 0.70; % this defines the intensity of grey. 1 is white; 0 is black.
grey = [c c c];

donLen = length(dat);
if len > donLen
    len = donLen;
end;

drng=[start:start+len-1];
%bar(data(rng,2),'histc');

trng=(0:len-1)*tstep;
ymax = max(dat(drng))*fac;

% plot the donor
bar1=bar(trng,dat(drng),'histc');
set(bar1,'FaceColor', 'r', 'EdgeColor', 'r');

%pbaspect([2.5 1 1]);

%set(gca,'XTick',-pi:pi/2:pi
%axis tight;
%xlim([1 len-1]);
%ylim([-ymax*1.2 ymax*1.2]);
set(gca,'FontName','Helvetica','FontSize',14);
ylim([0 ymax]);
xlim([0 trng(end)]);
ylabel(str,'FontName','Helvetica','FontSize',14);
xlabel('Time [ms]');
