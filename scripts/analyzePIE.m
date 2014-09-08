% FILENAME: analyzePIE.m, First Created on Jan. 04, 2012
%
% DESCRIPTION: This is the master file for reading in and processing a .pt3
% file
%
% DEPENDENT FILES: 
%   readHeader
%   readCounts2
%
%
% INPUT FILES15

% 
%
% OUTPUT FILES


% INSTRUCTIONS:
%   Set USER VARIABLES (filename, pathname) either using the gui or,
%   preferably, by defining them explicitly. Your choice is registered
%   using GUI_PREF

% BACKGROUND INFORMATION
%
%   Analyze the buffer first. This gives us the time gate thresholds (from 
%   TCSPC) and an estimate of background counts (from MCS).


% makeMCS.m
% This script reads a PicoHarp 300 T3 Mode data file (*.pt3)
% Works with file format version 2.0 only!

% Tested with Matlab R2011b
% John Robinson, SDSU, updated January 2012
% this script was derived from read_pt3.m 
% Peter Kapusta, Michael Wahl, PicoQuant GmbH 2007, updated May 2007

% NOTES
% Note that marker events have a lower time resolution and may therefore appear 
% in the file slightly out of order with respect to regular (photon) event records.
% This is by design. Markers are designed only for relatively coarse 
% synchronization requirements such as image scanning. 

% JMR: we don't care about markers

% T3 Mode data are written to an output file [filename].out 
% We do not keep it in memory because of the huge amout of memory
% this would take in case of large files. Of course you can change this, 
% e.g. if your files are not too big. 
% Otherwise it is best process the data on the fly and keep only the results.

% JMR: I have substantially rewritten read_pt3.m as the function readCounts2.
% Time speed up is considerable. Memory is not an issue. 

% BEFORE WE BEGIN

clear all;
close all;
clc;

% BEGIN

% USER PREFERENCES 

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

PLOT_VS_TIME = 0 % plot TCSPC vs (1) time & (0) channel #
GUI_PREF = 1; % 1 = true; 0 = false
COMMENT = 1; % add comments to graph. Supress for publication.

%PLOTS
PLOT_PHOTON_COUNTING_HISTOGRAM = 0; % render the photon counting histograms
PLOT_BURST = 0;
PLOT_TCSPC = 0; %
PLOT_BURST_SEL = 0;
PLOT_1D_HIST = 1; % 

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

% Set filename and pathname using the GUI

if (GUI_PREF == 1)
    [filename, pathname]=uigetfile('*.pt3', 'T3 Mode data:');
else
% Set filename and pathname manually
    %filename = 'D11A647(50pM PIE).pt3';
    %pathname = '../input/';
end;
inFile=fopen([pathname filename]);
logFileName = [pathname filename(1:length(filename)-4) '.log'];
logFile = fopen(logFileName,'W');
fprintf(logFile,'Header from file: %s\n', filename);
fprintf(logFile,'\n');
[numRecords resolution syncPeriod] = readHeader(inFile, logFile);
[Ddt Dtt Adt Att] = readCounts2(inFile, logFile, numRecords, resolution, syncPeriod);
%[Ddt Dtt Adt Att] = readCounts2(inFile, logFile, 10000000, resolution, syncPeriod);
fclose(inFile);
fclose(logFile);

%%
% ANALYSIS-DEPENDENT USER VARIABLES 
% 
% ------------------- User selectable variables ------------
% 
% %%
% % gamma = input('gamma value:  ');
PIE_MIN = input('PIE minimum value(~15):  ');
% PIE_MAX = input('pie maximum value:  ');
NOISE_MIN = input('noise minimum value(~10):  ');
% NOISE_MAX = input('noise maximum value:  ');

INCLUDE_SCATTERING = 1; % this shifts the time gate back to include the whole TCSPC histogram
BIN_TIME = 1.0; % in millisecond
% PIE_MIN = 15;   % input value
PIE_MAX = 70;
% NOISE_MIN = 10;  % input value
NOISE_MAX = 70;
%GAMMA = 2.1; %for Troponin (Cy3-Atto655)
GAMMA = 1.15; % for DNA (Cy3-Alexa647)

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%DESC = 'scattering included';
%DESC = 'D11A647(50pM), 1 msec';
%NAME = input('Description of output data:  ');
%DESC = 'D7A647(50pM)'; %$$$$$$$$$$$$$$$ Type data file name $$$$$$$$$$$$$$
%filename = filename(length(filename)-4:length(filename)); %GHK
%
filename = filename(1:length(filename)-4); %GHK for printing file
DESC = filename;                           %GHK
if not(length(DESC))
    DESC_LABEL = '';
else
    DESC_LABEL = sprintf('(%s) ',DESC);
end;
%
%% =========================== TIME GATING ===============================
%  forget about trying to get rid of scattered light.
%
load '../input/timeGateFilters.mat' DDTimeGate DATimeGate AATimeGate; % these are from analyzeTimeGates.m
%
FWHM = 7;
if INCLUDE_SCATTERING
    DDTimeGate(1) = floor(DDTimeGate(1) - 2.5*FWHM);
    DATimeGate(1) = floor(DATimeGate(1) - 2.5*FWHM);
    AATimeGate(1) = floor(AATimeGate(1) - 2.5*FWHM);
end;
%
% $$$$$$ Time gate set by manual, GHK
%
%DDTimeGate = [100, 400]; %[100, 800]
%DATimeGate = [115, 500]; %[115, 800]
%AATimeGate = [900, 1300]; %[920, 1500]
%
% TCSPC plots in terms of time
%
chan = 0:round(syncPeriod/resolution); % gives integers
%
donorHist = histc(Ddt,chan); % do it on the center instead of edges (hist vs. histc)
acceptorHist = histc(Adt,chan);
% USER VARIABLES
%
spadShift = DDTimeGate(1) - DATimeGate(1); % just for fun

if PLOT_VS_TIME
    t = chan*resolution; % in nsec
    % green excitation
    subplot(2,1,1);
    rectangle('Position',[resolution*DDTimeGate(1),0,resolution*(DDTimeGate(2)-DDTimeGate(1)),max(donorHist)*1.1],'FaceColor',[0.9 0.9 0.9],'edgecolor','none');
    hold on;
    plot(t,donorHist,'b');
    % semilogy(timeChan,donorHist,'b');
    hold off;
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',14);
    ylabel('Counts');
    xlabel('Time [nsec]'); 
    title('Donor');
    grid on;

    subplot(2,1,2);
    rectangle('Position',[resolution*DATimeGate(1),0,resolution*(DATimeGate(2)-DATimeGate(1)),max(acceptorHist)*1.1],'FaceColor',[0.9 0.9 0.9],'edgecolor','none');
    hold on;
    rectangle('Position',[resolution*AATimeGate(1),0,resolution*(AATimeGate(2)-AATimeGate(1)),max(acceptorHist)*1.1],'FaceColor',[0.9 0.9 0.9],'edgecolor','none');
    plot(t,acceptorHist,'r');
    hold off;
    % semilogy(timeChan,acceptorHist,'r');
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',14);
    ylabel('Counts'); 
    xlabel('Time [nsec]'); 
    title('Acceptor');
    grid on;

    str = sprintf('../output/%sTCSPC_Time.eps', DESC_LABEL); print('-depsc',str)
    if PLOT_TCSPC
        figure;
    else
        clf;
    end;

else % plot vs. channel
            % green excitation
    subplot(2,1,1);
    rectangle('Position',[DDTimeGate(1),0,(DDTimeGate(2)-DDTimeGate(1)),max(donorHist)*1.1],'FaceColor',[0.9 0.9 0.9],'edgecolor','none');
    hold on;
    plot(chan,donorHist,'b');
    % semilogy(timeChan,donorHist,'b');
    hold off;
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',14);
    ylabel('Counts');
    xlabel('Channel'); 
    title('Donor');
    grid on;

    subplot(2,1,2);
    rectangle('Position',[DATimeGate(1),0,(DATimeGate(2)-DATimeGate(1)),max(acceptorHist)*1.1],'FaceColor',[0.9 0.9 0.9],'edgecolor','none');
    hold on;
    rectangle('Position',[AATimeGate(1),0,(AATimeGate(2)-AATimeGate(1)),max(acceptorHist)*1.1],'FaceColor',[0.9 0.9 0.9],'edgecolor','none');
    plot(chan,acceptorHist,'r');
    hold off;
    % semilogy(timeChan,acceptorHist,'r');
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',14);
    ylabel('Counts');
    xlabel('Channel'); 
    title('Acceptor');
    grid on;
    %print -depsc2 '../output/TCSPCchan.eps'
    str = sprintf('../output/%sTCSPC_chan.eps', DESC_LABEL); print('-depsc',str)
    if PLOT_TCSPC
       figure;
    else
        clf;
    end;
end;

%%
% ============================ MCS =====================================

% MCS Plots with and without the time gate.

BIN_TIME_INC = BIN_TIME * 1E6; % this gives it to us in nsec.
% Apply time gates
DDFilter = (Ddt >= DDTimeGate(1)) & (Ddt <= DDTimeGate(2)); 
DAFilter = (Adt >= DATimeGate(1)) & (Adt <= DATimeGate(2));
AAFilter = (Adt >= AATimeGate(1)) & (Adt <= AATimeGate(2));

maxTime = min([Dtt(end);Att(end)]/1E6); % in milliseconds % divide by 60,000 to get minutes

DDtmp = getMCS(Dtt,DDFilter, BIN_TIME,maxTime);
DAtmp = getMCS(Att,DAFilter, BIN_TIME,maxTime);
AAtmp = getMCS(Att,AAFilter, BIN_TIME,maxTime);

% note that length(DD) = length(DA) = length(AA)

smallest = min([length(DDtmp), length(DAtmp), length(AAtmp)]);

DD = DDtmp(1:smallest);
DA = DAtmp(1:smallest);
AA = AAtmp(1:smallest);

meanCounts = [mean(DD), mean(DA), mean(AA)]; % gives average counts per bin
stdCounts = [std(DD), std(DA), std(AA)];

if PLOT_BURST
    figure;
    start = 1000;
    len = 2000;
    str = sprintf('Counts / %1.1f ms',BIN_TIME)
    plotBurst1(BIN_TIME,DD,start,len,str); 
    title('DD');
    figure;
    plotBurst1(BIN_TIME,DA,start,len,str);
    title('DA');
    figure;
    plotBurst1(BIN_TIME,AA,start,len,str);
    title('AA');
end;

% Plots of the photon counting histograms. This shows the effect of the thresholds
% ====================== DD photon counting histogram ============
countBins = 0:max(DD);
DDCounts = myHistc(DD,countBins)/maxTime;
semilogx(DDCounts,countBins,'b');
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
xlabel('Frequency [Events/sec]');
ylab = sprintf('Counts / %1.1f ms',BIN_TIME); ylabel(ylab); 
title('DD photon counting histogram');
grid on;
ax = axis;
hold on;
cMin = 1E-8;
cMax = max(DDCounts);
semilogx([cMin cMax], [NOISE_MIN NOISE_MIN],'r');
semilogx([cMin cMax], [NOISE_MAX NOISE_MAX],'g');
%xlim([cMin cMax]);
hold off;
axis(ax);
% print -depsc2 '../output/DDphotonCountingHistogram.eps'
str = sprintf('../output/%sDDphotonCountingHistogram.eps', DESC_LABEL); print('-depsc',str);
if PLOT_PHOTON_COUNTING_HISTOGRAM
    figure;
else
    clf;
end;

% ====================== DA photon counting histogram ============
countBins = 0:max(DA);
DACounts = myHistc(DA,countBins)/maxTime;
semilogx(DACounts,countBins,'b');
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
xlabel('Frequency [Events/sec]');
ylab = sprintf('Counts / %1.1f ms',BIN_TIME); ylabel(ylab); 
title('DA photon counting histogram');
grid on;
% ax = axis;
% hold on;
% cMin = 1E-8;
% cMax = max(DACounts);
% semilogx([cMin cMax], [NOISE_MIN NOISE_MIN],'r');
% semilogx([cMin cMax], [NOISE_MAX NOISE_MAX],'g');
% %xlim([cMin cMax]);
% hold off;
% axis(ax);
% print -depsc2 '../output/DDphotonCountingHistogram.eps'
str = sprintf('../output/%sDAphotonCountingHistogram.eps', DESC_LABEL); print('-depsc',str);
if PLOT_PHOTON_COUNTING_HISTOGRAM
    figure;
else
    clf;
end;

% ====================== DD + DA photon counting histogram ============
DTot = DD + DA;
countBins = 0:max(DTot);
DDDACounts = myHistc(DTot,countBins)/maxTime;
% figure
semilogx(DDDACounts,countBins,'b');
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
xlabel('Frequency [Events/sec]');
ylabel(ylab);
title('DD + DA photon counting histogram');
grid on;
ax = axis;

hold on;
%myLim = xlim;
cMin = 1E-8;
%cMin = min(DDDAcounts)
cMax = max(DDDACounts);
semilogx([cMin cMax], [NOISE_MIN NOISE_MIN],'r');
semilogx([cMin cMax], [NOISE_MAX NOISE_MAX],'g');
% xlim([cMin cMax]);
hold off;
%axis tight;
axis(ax);
% print -depsc2 '../output/DDDAphotonCountingHistogram.eps'
str = sprintf('../output/%sDDDAphotonCountingHistogram.eps', DESC_LABEL); print('-depsc',str);
if PLOT_PHOTON_COUNTING_HISTOGRAM
    figure;
else
    clf;
end;

% ====================== AA photon counting histogram ============
countBins = 0:max(AA);
AACounts = myHistc(AA,countBins)/maxTime;
%bar(x,log10(n));
% figure
semilogx(AACounts,countBins,'b');
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
xlabel('Frequency [Events/sec]');
ylabel(ylab);
title('AA photon counting histogram');
grid on;
ax = axis;
hold on;
cMin = 1E-8;
cMax = max(AACounts);
semilogx([cMin cMax], [PIE_MIN PIE_MIN],'r');
semilogx([cMin cMax], [PIE_MAX PIE_MAX],'g');
%xlim([cMin cMax]);
hold off;
axis(ax);
% print -depsc2 '../output/AAphotonCountingHistogram.eps'
str = sprintf('../output/%sAAphotonCountingHistogram.eps', DESC_LABEL); print('-depsc',str);
if PLOT_PHOTON_COUNTING_HISTOGRAM
    figure;
else
    clf;
end;

% ============================= Main analysis =============================

NUM_BINS = 41;
TE_LIM = [0.0 1.0];
S_LIM = [0.0 1.0];
TEBins = linspace(TE_LIM(1),TE_LIM(2),NUM_BINS);
SBins = linspace(S_LIM(1),S_LIM(2),NUM_BINS);
%bins = (0:numBins)/double(numBins);
%bins = (0:(numBins-1))/(numBins-1.);

start = 200000;
len = 3000;
len = len + 1; %makes the plots prettier.

% =============== Calculate the filters ======================
% define PIE filter
redFilter = ((AA > PIE_MIN) & (AA < PIE_MAX)); %pieFilter is a logical array
greenFilter = ((DD+DA > NOISE_MIN) & (DD+DA < NOISE_MAX));
%countFilter = ((DD > noiseMin) & (DD < noiseMax));
greenAndRedFilter = redFilter & greenFilter;

greenOrRedFilter = redFilter | greenFilter;


% ------------------- Burst Traces ----------------------------------

% ------------------- without PIE filtering ------------------------
%figure
plotBurst2sel(BIN_TIME,DD,DA,greenFilter,start,len);
pbaspect([4.0 1 1]);
title('Donor excitation'); %GHK
%print -depsc2 '../output/normalBurst.eps';
str = sprintf('../output/%sgreenExcitedBurst.eps', DESC_LABEL); print('-depsc',str);
if PLOT_BURST_SEL
    figure;
else
    clf;
end;

% ------------------- with PIE filtering ------------------------
% Plot results (as subplots)
% top is burst plot of donor and acceptor with Donor excitation
% figure
subplot(2,1,1);
%plotBurst2sel(tstep,DD,DA,pieFilter,start,len);
plotBurst2sel(BIN_TIME,DD,DA,greenAndRedFilter,start,len);
title('Donor excitation');

% set(h,'XTick',[])
% set(gca,'XTickLabel',[]);

% bottom plot is plot of acceptor getting directly excited
subplot(2,1,2); 
plotBurst1sel(BIN_TIME,AA,redFilter,start,len);
title('Acceptor excitation');
% pbaspect([4.0 1 1]);

str = sprintf('../output/%spieBurst.eps', DESC_LABEL); print('-depsc',str);
if PLOT_BURST_SEL
    figure;
else
    clf;
end;

%
%
%% ----------------------- 1D histograms ---------------------------------
% the 1D analysis with the PIE selection on.

F_D1 = GAMMA*DD + DA;
F_D2 = DD + DA; 
F_A = AA;

TE = DA./F_D1;
S = F_D2./(F_D2 + F_A); % S used a gamma uncorrected DD.

% fTE = filterSelect(TE,greenFilter); this selects the green-excited
fTE = filterSelect(TE,greenOrRedFilter);
fTEPIE = filterSelect(TE,greenAndRedFilter); %with PIE selection
TEHist = myHistc(fTE,TEBins);
TEHistPIE = myHistc(fTEPIE,TEBins); %with PIE selection

% fS = filterSelect(S,greenFilter);
fS = filterSelect(S,greenOrRedFilter);
fSPIE = filterSelect(S,greenAndRedFilter);
SHist = myHistc(fS,SBins);
SHistPIE = myHistc(fSPIE,SBins);

% ------- Plot of TE
figure(1),
subplot(2,1,1);
TEPlotPIE = bar(TEBins,TEHistPIE,1,'histc');
set(TEPlotPIE,'EdgeColor','k');
set(TEPlotPIE,'Facecolor','w');
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
%xlim([TEBins(1) TEBins(end)]);
xlim(TE_LIM);
xlabel('TE'); ylabel('Counts');
title('With PIE selection');

subplot(2,1,2);
TEPlot = bar(TEBins,TEHist,1,'histc');
set(TEPlot,'EdgeColor','k');
set(TEPlot,'Facecolor','w');
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
%xlim([TEBins(1) TEBins(end)]);
xlim(TE_LIM);
xlabel('TE'); ylabel('Counts');
title('Without PIE selection');
if COMMENT
    legend(DESC);
end;
str = sprintf('../output/%s1D_TE.eps', DESC_LABEL); print('-depsc',str);
if PLOT_1D_HIST
    figure;
else
    clf;
end;

% ----- Plot of S -------
figure(2),
subplot(2,1,1);
SPlotPIE = bar(SBins,SHistPIE,1,'histc');
set(SPlotPIE,'EdgeColor','k');
set(SPlotPIE,'Facecolor','w');
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
%xlim([SBins(1) SBins(end)]);
xlim(S_LIM);
xlabel('S'); ylabel('Counts'); 
title('With PIE selection');

subplot(2,1,2);
SPlot = bar(SBins,SHist,1,'histc');
set(SPlot,'EdgeColor','k');
set(SPlot,'Facecolor','w');
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
%xlim([SBins(1) SBins(end)]);
xlim(S_LIM);
xlabel('S'); ylabel('Counts'); 
title('Without PIE selection');
if COMMENT
    legend(DESC);
end;
str = sprintf('../output/%s1D_S.eps', DESC_LABEL); print('-depsc',str);
if PLOT_1D_HIST
    figure;
else
    clf;
end;


%% -------------------------2D Hist Log scale in 2D hist-------------------
% for 2D histogram analysis, we select only based on countFilter
figure(3),
subplot(2,2,1); 
bar(TEBins,TEHist,1,'histc'); 
TEHistPlot = gca; 
axis([TE_LIM 0 max(TEHist)*1.01]); %axis('off');
%set(TEHistPlot,'xtick',[]); % don't write the x axis
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
ylabel('Counts');
title('2D histogram (Log) w/o PIE selection');
% grid on;
% Main plot
subplot(2,2,3); 
%mainHist(fTE,fS,TEBins,SBins,1,'k.',flipud(gray(256)));
% mainHistLog(fTE,fS,TEBins,SBins,1,'k.',flipud(bone(256)));  %JMR
mainHistLog(fTE,fS,TEBins,SBins,1,'k.',jet(256));  %% Log scale in 2D
h1 = gca; % axis([xlim ylim]);
xlabel('TE'); ylabel('S');

% SHistPlot
subplot(2,2,4);     
barh(SBins,SHist,1,'histc');
SHistPlot = gca; 
axis([0 max(SHist)*1.01 S_LIM]); %axis('off');
%set(SHistPlot,'ytick',[]); % don't write the y axis
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
xlabel('Counts');
% grid on;
%line([0 0],ylim-yoff,'Color','k')

%set(h1,'Position',[0.1 0.1 0.60 0.60]);
%set(TEHistPlot,'Position',[.1 .75 .60 .2]);
%set(SHistPlot,'Position',[.75 .1 .2 .60]);
set(h1,'Position',[0.15 0.1 0.50 0.60]);
set(TEHistPlot,'Position',[.15 .75 .50 .2]);
set(SHistPlot,'Position',[.70 .1 .2 .60]);
str = sprintf('../output/%s2DHistThreshLog.eps', DESC_LABEL); print('-depsc',str);

if COMMENT
    aa = legend(DESC);
    set (aa, 'Position',[.70 .82 .20 .05]);
end;
%str = sprintf('../output/%s2DHistThresh.eps', DESC_LABEL); print('-depsc',str);


%% -------------------------2D Hist ------------------------------------
% for 2D histogram analysis, we select only based on countFilter
figure(4),
subplot(2,2,1); 
bar(TEBins,TEHist,1,'histc'); 
TEHistPlot = gca; 
axis([TE_LIM 0 max(TEHist)*1.01]); %axis('off');
%set(TEHistPlot,'xtick',[]); % don't write the x axis
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
ylabel('Counts');
title('2D histogram w/o PIE selection');
% grid on;
% Main plot
subplot(2,2,3); 
%mainHist(fTE,fS,TEBins,SBins,1,'k.',jet(256));
%mainHist(fTE,fS,TEBins,SBins,1,'k.',flipud(bone(256)));
mainHist(fTE,fS,TEBins,SBins,1,'k.',flipud(hot(256)));  % Hot color
h1 = gca; % axis([xlim ylim]);
xlabel('TE'); ylabel('S');

% SHistPlot
subplot(2,2,4);     
barh(SBins,SHist,1,'histc');
SHistPlot = gca; 
axis([0 max(SHist)*1.01 S_LIM]); %axis('off');
%set(SHistPlot,'ytick',[]); % don't write the y axis
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
xlabel('Counts');
% grid on;
%line([0 0],ylim-yoff,'Color','k')

%set(h1,'Position',[0.1 0.1 0.60 0.60]);
%set(TEHistPlot,'Position',[.1 .75 .60 .2]);
%set(SHistPlot,'Position',[.75 .1 .2 .60]);
set(h1,'Position',[0.15 0.1 0.50 0.60]);
set(TEHistPlot,'Position',[.15 .75 .50 .2]);
set(SHistPlot,'Position',[.70 .1 .2 .60]);
str = sprintf('../output/%s2DHistThresh.eps', DESC_LABEL); print('-depsc',str);

if COMMENT
    bb = legend(DESC);
    set (bb, 'Position',[.70 .82 .20 .05]);
end;
%str = sprintf('../output/%s2DHistThresh.eps', DESC_LABEL); print('-depsc',str);


%% for 2D histogram with PIE filtering %GHK

figure(5),
subplot(2,2,1); 
bar(TEBins,TEHistPIE,1,'histc'); 
TEHistPlot = gca; 
axis([TE_LIM 0 max(TEHistPIE)*1.01]); %axis('off');
%set(TEHistPlot,'xtick',[]); % don't write the x axis
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
ylabel('Counts');
title('2D histogram with PIE selection');

% grid on;
% Main plot
subplot(2,2,3); 
%mainHist(fTEPIE,fSPIE,TEBins,SBins,1,'k.',flipud(gray(256)));
%mainHist(fTEPIE,fSPIE,TEBins,SBins,1,'k.',flipud(bone(256)));  %JMR
mainHist(fTEPIE,fSPIE,TEBins,SBins,1,'k.',jet(256));  %%MMK
h1 = gca; % axis([xlim ylim]);
xlabel('TE'); ylabel('S');

% SHistPlot
subplot(2,2,4);     
barh(SBins,SHistPIE,1,'histc');
SHistPlot = gca; 
axis([0 max(SHistPIE)*1.01 S_LIM]); %axis('off');
%set(SHistPlot,'ytick',[]); % don't write the y axis
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
xlabel('Counts');
% grid on;
%line([0 0],ylim-yoff,'Color','k')

%set(h1,'Position',[0.1 0.1 0.60 0.60]);
%set(TEHistPlot,'Position',[.1 .75 .60 .2]);
%set(SHistPlot,'Position',[.75 .1 .2 .60]);
set(h1,'Position',[0.15 0.1 0.50 0.60]);
set(TEHistPlot,'Position',[.15 .75 .50 .2]);
set(SHistPlot,'Position',[.70 .1 .2 .60]);

str = sprintf('../output/%s2DHistThreshPIE.eps', DESC_LABEL); print('-depsc',str);

if COMMENT
   cc = legend(DESC);
   set (cc, 'Position',[.70 .82 .20 .05]);
end;

%str = sprintf('../output/%s2DHistThreshPIE.eps', DESC_LABEL); print('-depsc',str);

%% for peak point of TE, %GHK, Jan. 29, 2012

figure(6),
subplot(2,1,1)
%scs = csaps (TEBins, TEHistPIE);
%fnplt (scs, 2, 'r'), 
xxx=0:0.025:1;
yyy=spline(TEBins, TEHistPIE, xxx);
%plot (TEBins, TEHistPIE,'or',xxx,yyy,'r');
plot (xxx,yyy,'r','linewidth',2);
hold on;
bar(TEBins, TEHistPIE, 'b');
hold off;
%
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
xlim(TE_LIM);
xlabel('TE'); ylabel('Counts');
title(['TE of ',filename]);
%title('TE With PIE selection');
%if COMMENT
%   legend(DESC);
%end;
%
% cubic smoothing spline of TE
subplot(2,1,2)
scs = csaps (TEBins, TEHistPIE);
fnplt (scs, 2, 'r'),
xlabel('TE');
grid on;
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
title('TE after smoothing');
%
%
%% use matlab curvefitting toolbox to fit to gaussians, Feb. 09 2012
%
% Here, what we are doing is doing histogrms on center instead of on edges.
%
% PARAMETERS
figure(7),
NUM_GAUSSIANS = 1;     % 2 for two gaussian fit
REDUCED_FIT_RANGE = 1; % if true the don't fit to the full range of TE;s
%                      % 0 for full range, 1 for reduced fit range 
%
R_NUM_BINS = NUM_BINS - 1; % doing on centers reduces NUM_BINS by 1
binWidthD2 = (TEBins(2)-TEBins(1))/2;
TEMidBins = linspace(TE_LIM(1)+binWidthD2,TE_LIM(2)-binWidthD2,R_NUM_BINS);
%
if REDUCED_FIT_RANGE
%TEFitRange = [0.5 0.98]; % reduced fit range for 0.8 above
TEFitRange = [0.2 0.8]; % center weighted fitting
%TEFitRange = [0.02 0.55]; % reduced fit range for 0.3 below
else
   TEFitRange = TE_LIM;
%   TEFitRange = [0.1 0.9];
end
% TEFitRange = [0 1];
% alternative is to use function excludedata
lowerIndex = min(find(TEMidBins > TEFitRange(1)));
upperIndex = max(find(TEMidBins < TEFitRange(2)));
% So, now find the indeces in TEBins(2) for the lims of TERange.
%
% plot just the data in the fitting range
% bar(TEMidBins(lowerIndex:upperIndex)',TEHistPIE(lowerIndex:upperIndex));
%
if (NUM_GAUSSIANS > 3) 
    return;
end;

switch NUM_GAUSSIANS
    case 1
        ft = fittype('gauss1');
    case 2
        ft = fittype('gauss2');
    case 3
        ft = fittype('gauss3');
end
%
[a,b,c] = fit(TEMidBins(lowerIndex:upperIndex)',TEHistPIE(lowerIndex:upperIndex),ft); %
%[a,b,c] = fit(TEMidBins(lowerIndex:upperIndex)',TEHistPIE(lowerIndex:upperIndex),ft,'startpoint',[100,1.0,0.1]); %
methods(a);
bar(TEMidBins,TEHistPIE(1:R_NUM_BINS));
hold on;
fp = plot(a,'-r');
hold off;
set(fp,'linewidth',2);
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
% xlim(TE_LIM);
xlabel('TE'); ylabel('Counts');
title(['TE of ',filename]);

% argnames(a);
c = coeffvalues(a); 
meanTE = c(2);
stdTE = c(3);
%
switch NUM_GAUSSIANS
    case 1
        meanTE = c(2);
    case 2
        meanTE = (c(1)*c(2) + c(4)*c(5) )/( c(1) + c(4) ); % amplitude weighted mixture
    case 3
        meanTE = (c(1)*c(2) + c(4)*c(5) + c(7)*c(8) )/( c(1) + c(4) + c(7) ); % amplitude weighted mixture
end
PIE_MIN, NOISE_MIN,
meanTE
MTE = meanTE; % for texting mean TE on figure
ytext = max(TEHistPIE);
text(0.01,ytext,['Mean TE = ', num2str(MTE)],'FontSize',14,'BackgroundColor',[.7 .9 .7]);
% for print output folder
str = sprintf('../output/%sGaussianFittedMeanTE.eps', DESC_LABEL); print('-depsc',str); %print on output folder
%
%
%
%
%%
%
figure(8),
%
bar(TEMidBins,TEHistPIE(1:R_NUM_BINS));
hold on;
fp = plot(a,'-r');
hold off;
set(fp,'linewidth',2);
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
% xlim(TE_LIM);
xlabel('TE'); ylabel('Counts');
title(['TE of ',filename]);

% argnames(a);
c = coeffvalues(a); 
meanTE = c(2);
stdTE = c(3);
Totalevents = sum(TEHistPIE); % sum of all FRET events after threshold
Totalevents,
switch NUM_GAUSSIANS
    case 1
        meanTE = c(2);
    case 2
        meanTE = (c(1)*c(2) + c(4)*c(5) )/( c(1) + c(4) ); % amplitude weighted mixture
    case 3
        meanTE = (c(1)*c(2) + c(4)*c(5) + c(7)*c(8) )/( c(1) + c(4) + c(7) ); % amplitude weighted mixture
end
%gtext(['TE = ', num2str(MTE)],'FontSize',14,'BackgroundColor',[.7 .9 .7]);
ytext = max(TEHistPIE);
yytext = ytext*0.1;
yyytext = ytext*0.9;
text(0.01,ytext,['Mean TE = ', num2str(MTE)],'FontSize',14,'BackgroundColor',[.7 .9 .7]);
text(0.4,yytext,['FRET events = ', num2str(Totalevents)],'FontSize',14,'BackgroundColor',[.7 .7 .7]);
text(0.01,yyytext, ['PIE Min: ', num2str(PIE_MIN), ',  Noise Min: ', num2str(NOISE_MIN)],'FontSize',10,'BackgroundColor',[.7 .9 .7]);
%gtext(['PIE Min: ', num2str(PIE_MIN), ', Noise Min: ', num2str(NOISE_MIN)],'FontSize',12,'BackgroundColor',[.7 .9 .7]);
grid
%
str = sprintf('../output/%sGaussianFittedMeanTE2.eps', DESC_LABEL); print('-depsc',str); %print on output folder
%
%
%