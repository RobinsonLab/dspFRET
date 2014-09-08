function plotBurst2sel(tstep,don,acc,sel,PP)
% plots both donor and acceptor channels in a burst trace. 
% those bursts that are selected by the PIE filter are shown with a grey
% background.

% INPUTS:
% tstep, binning time step [msec]
% don, a array of size (t,2) 
% acc, a array of size (t,2) 
% sel, selection filter, logical of length t
% PP, PLOT_PREFS (provides START and LEN)


% OUTPUTS:



% data is an t,2 array.
% sel is the select filter
% figure;

% CONSTANTS
fac = 1.2;
c = 0.70; % this defines the intensity of grey. 1 is white; 0 is black.
grey = [c c c];

drng = [PP.START:(PP.START + PP.LEN -1)];
%bar(data(rng,2),'histc');

trng=(0:(PP.LEN-1))*tstep;

maxval1 = max(don(drng));
maxval2 = max(acc(drng));
ymax = max([maxval1 maxval2])*fac;

% color code the PIE-selected bursts
barb1=bar(trng,sel(drng)*ymax,'histc'); % will need to modify to sel(drng,2)*ymax
set(barb1,'FaceColor', grey, 'EdgeColor', grey);
hold on;
barb2=bar(trng,-sel(drng)*ymax,'histc'); % will need to modify to -sel(drng,2)*ymax
set(barb2,'FaceColor', grey, 'EdgeColor', grey);

% plot the donor
bar1=bar(trng,don(drng),'histc');
set(bar1,'FaceColor', 'b', 'EdgeColor', 'b');

% plot the acceptor
bar2=bar(trng,-acc(drng));
set(bar2,'FaceColor', 'r', 'EdgeColor', 'r');

%plot(trng,sel*ymax*1.1,'.k');

hold off;
pbaspect([4.0 1 1]);
%set(gca,'XTick',-pi:pi/2:pi
%axis tight;
%xlim([1 len-1]);
%ylim([-ymax*1.2 ymax*1.2]);
ylim([-ymax ymax]);
xlim([0 trng(end)]);
set(gca,'FontName','Helvetica','FontSize',14);
ylabel('Counts');
xlabel('Time [msec]');
