function plotBurst1sel(tstep,don,sel,PP)
% plots both donor and acceptor channels in a burst trace. 
% those bursts that are selected by the PIE filter are shown with a grey
% background.

% INPUTS:
% tstep, binning time step [msec]
% don, a array of size (t,2) 
% sel, selection filter, logical of length t
% start, start channel to plot data
% len, number of bins to plot


% OUTPUTS:


% CONSTANTS
fac = 1.2;
c = 0.70; % this defines the intensity of grey. 1 is white; 0 is black.
grey = [c c c];

drng = [PP.START:(PP.START + PP.LEN -1)];
%drng=[start:start+len-1];
%bar(data(rng,2),'histc');

trng=(0:(PP.LEN-1))*tstep;
ymax = max(don(drng))*fac;


% color code the PIE-selected bursts
barb1=bar(trng,sel(drng)*ymax,'histc'); % will need to modify to sel(drng,2)*ymax
set(barb1,'FaceColor', grey, 'EdgeColor', grey);
hold on;
% plot the donor
bar1=bar(trng,don(drng),'histc');
set(bar1,'FaceColor', 'r', 'EdgeColor', 'r');
hold off;
%pbaspect([2.5 1 1]);

%set(gca,'XTick',-pi:pi/2:pi
%axis tight;
%xlim([1 len-1]);
%ylim([-ymax*1.2 ymax*1.2]);
set(gca,'FontName','Helvetica','FontSize',14);
ylim([0 ymax]);
xlim([0 trng(end)]);
ylabel('Counts','FontName','Helvetica','FontSize',14);
xlabel('Time [msec]');
