% FILENAME: analyzePIE.m, First Created on Jan. 04, 2012
%
% REVISION HISTORY:
%
% Ver  Author   Date    Description 
% ---  ------   ----    --------------------------------------------
%   1   JMR     1/4/12  o Read in TTRT file
%                       o PCH
%                       o 2DFRET hist.

%   2   JMR     6/8/12  o Lifetime anslysis of donor-only and donor-acceptor 
%                       o put variables into data structures
%                       o lots of comments
%                       o maps for photon sorting into DD, DA, AA 
%                       o maps for assigning photons to bins

% DESCRIPTION: 
% This is the master file for reading in and processing single molecuele FRET data
%
% DEPENDENT FILES: 
%   readHeader      reads in the header from a TTTR file
%   readCounts2     reads the contents of a TTTR file
%
%
% INPUT FILES:
% TTTR file in ../input
% 
%
% OUTPUT FILES:
%  ../output

% 1) figures: 
%
%
%
%
% 2) log file:
% ../output/filename.log
%
%

% USER INSTRUCTIONS:
%   Set PREFS and PARAMS

% BACKGROUND INFORMATION:
%
%   Analyze the buffer first. This gives us the time gate thresholds (from 
%   TCSPC) and an estimate of background counts (from MCS).



% BEFORE WE BEGIN

clear all;
close all;
clc;

% BEGIN

% USER PREFERENCES 
PREFS.GUI = 1; % 1 = true; 0 = false

% Render plots?
PREFS.PLOT_VS_TIME = 1; % plot TCSPC vs (1) time & (0) channel #
PREFS.COMMENT = 1; % add comments to graph. Supress for publication.
PREFS.PLOT_PHOTON_COUNTING_HISTOGRAM = 0; % render the photon counting histograms
PREFS.PLOT_BURST = 0;
PREFS.PLOT_TCSPC = 0; %
PREFS.PLOT_BURST_SEL = 0;
PREFS.PLOT_1D_HIST = 1; % 
PREFS.VERBOSE = 1; % write out what is going on.

% DEFAULT PARAMS
PARAMS.INCLUDE_SCATTERING = 0; % includes scattering that is evident as the spike in the TCSPC histogram?
PARAMS.BIN_TIME = 1.0; % in milliseconds
PIE_MIN = 15;
PARAMS.PIE_MAX = 70;
NOISE_MIN = 10;
PARAMS.NOISE_MAX = 70;
%GAMMA = 2.1; %for Troponin (Cy3-Atto655)
% PARAMS.GAMMA = 1.15; % for DNA (Cy3-Alexa647)
PARAMS.GAMMA = 1.0; 
PARAMS.NOT_RED = 3; % this is very tight. Used to identify a donor-only species. Want to make sure that acceptor is definately not there!
PARAMS.NOT_GREEN = 3; % this is very tight. Used to identify a acceptor-only species. Want to make sure that donor is definately not there!

% PLOTTING PREFS
PLOT_PREFS.NUM_BINS = 41;
PLOT_PREFS.TE_LIM = [0.0 1.0];
PLOT_PREFS.S_LIM = [0.0 1.0];
% start = 200000; % or whatever....
PLOT_PREFS.START = 10;
PLOT_PREFS.LEN = 3000;
PLOT_PREFS.LEN = (PLOT_PREFS.LEN + 1); %makes the plots prettier.

%
% ================= READ input and set params ======================
%

%
if (PREFS.GUI == 1)
% Set filename and pathname using the GUI
    [inFile.name, inFile.path]=uigetfile('*.pt3', 'T3 Mode data:');
else
% Set filename and pathname manually
    inFile.name = 'D11A647(50pM PIE).pt3';
    inFile.path = '../input/';
end;
inFile.base = inFile.name(1:length(inFile.name)-4);
logFile.name = [inFile.path inFile.base '.log']; % note, this is how to cat strings

inFile.ptr = fopen([inFile.path inFile.name]);

% >> write out the header
logFile.ptr = fopen(logFile.name,'W');
fprintf(logFile.ptr,'\n---------------------------------------------------------\n');
fprintf(logFile.ptr,'>> Header from input file: %s\n', inFile.name);
fprintf(logFile.ptr,'\n---------------------------------------------------------\n');
[inFile.numRecords inFile.resolution inFile.syncPeriod] = readHeader(inFile.ptr, logFile.ptr);
fclose(logFile.ptr);

str = logFile.name
diary(str);

% >> get the photon events
% inFile.numRecords = 100000; % use this when debugging
[Dlc Dgt Alc Agt] = readCounts2(PREFS, inFile);
fclose(inFile.ptr);
% diary off;

% Main data here......
% For each count that is registered in the donor channel (spad), we have a local
% time of that count (Dlc) and the global time of that count (Dgt). Dlc is
% used for TCSPC and Dgt is used for MCS.

% For each count that is registered in the acceptor channel (spad), we have a local
% time of that count (Alc) and the global time of that count (Agt). Alc is
% used for TCSPC and Agt is used for MCS. We need to further process Alc to
% figure out whether the photon came from excitation with "green" laser (DA) or
% from the "red" laser (AA). This processing will be done later.

% Dlc, donor APD, local channel #
% Dgt, donor APD, global time
% Alc, acceptor APD, local channel #
% Agt, acceptor APD, global time

% Why are we reading into arrays?
% Would make sense to package these into an event structure
% Event
%  |
%  |-Donor- local- channel
%  |   |      \--- time
%  |   \--- global time
%  |
%  |-Acceptor- local- channel
%  |   |         \--- time
%  |   \--- global time
% 
% Then have an array of events

% But, we're going for speed,and accessing fields takes time?
 
% 
% ------------------- User selectable variables ------------
% 

fprintf(1,'---------------------------------------------------------\n');
fprintf(1,'>> Parameter selection\n');
fprintf(1,'---------------------------------------------------------\n');


% Override PIE_MIN and NOISE_MIN
if (PREFS.GUI == 1)
    fprintf(1, '\n');
    fprintf(1, 'Set thresholds');
    fprintf(1, '\n');
    pm = input(['AA minimum value(' int2str(PARAMS.PIE_MIN) ')']);
    if (length(pm) == 1)
        PARAMS.PIE_MIN = pm; % user didn't just hit enter
    end;
    clear pm;
    nm = input(['DA + AA minimum value(' int2str(PARAMS.NOISE_MIN) ')']);
    if (length(nm) == 1)
        PARAMS.NOISE_MIN = nm; % user didn't just hit enter
    end;
    clear nm;
end;

% This is a description that you may want to include in a plot (or other
% ouputs.

desc_str = inFile.base; % default is to use the pathname.                  
% other possibilities for DESC
if not(length(desc_str))
    PREFS.DESC_LABEL = '';
else
    PREFS.DESC_LABEL = sprintf('(%s) ',desc_str);
end;
%



% display(inFile)
fprintf(1,'---------------------------------------------------------\n');
fprintf(1,'>> Entering analysis stage with these parameters\n');
fprintf(1,'---------------------------------------------------------\n');

display(PREFS)
display(PARAMS)

%
% =================================================================
%   TCSPC Histogram
% =================================================================
% The idea was to try to get rid of
% scattered light by chopping off the "first part" of the TCSPC trace.

% This idea didn't really work out very well. 
%
load 'timeGateFilters.mat' DDTimeGate DATimeGate AATimeGate; % these are from analyzeTimeGates.m

% If you think that time gating is a bad idea and you want to include
% scattering then set PARAMS.INCLUDE_SCATTERING == 1

FWHM = 7;
if (PARAMS.INCLUDE_SCATTERING == 1)
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
chan = 0:round(inFile.syncPeriod/inFile.resolution); % gives integers
%
donorHist = histc(Dlc,chan); % do it on the center instead of edges (hist vs. histc)
acceptorHist = histc(Alc,chan);
% USER VARIABLES
%
spadShift = DDTimeGate(1) - DATimeGate(1); % just for fun

if (PREFS.PLOT_VS_TIME == 1)
    t = chan*inFile.resolution; % in nsec
    % green excitation
    subplot(2,1,1);
    rectangle('Position',[inFile.resolution*DDTimeGate(1),0,inFile.resolution*(DDTimeGate(2)-DDTimeGate(1)),max(donorHist)*1.1],'FaceColor',[0.9 0.9 0.9],'edgecolor','none');
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
    rectangle('Position',[inFile.resolution*DATimeGate(1),0,inFile.resolution*(DATimeGate(2)-DATimeGate(1)),max(acceptorHist)*1.1],'FaceColor',[0.9 0.9 0.9],'edgecolor','none');
    hold on;
    rectangle('Position',[inFile.resolution*AATimeGate(1),0,inFile.resolution*(AATimeGate(2)-AATimeGate(1)),max(acceptorHist)*1.1],'FaceColor',[0.9 0.9 0.9],'edgecolor','none');
    plot(t,acceptorHist,'r');
    hold off;
    % semilogy(timeChan,acceptorHist,'r');
    set(gca,'FontName','Helvetica');
    set(gca,'FontSize',14);
    ylabel('Counts'); 
    xlabel('Time [nsec]'); 
    title('Acceptor');
    grid on;

    str = sprintf('../output/%sTCSPC_Time.eps', PREFS.DESC_LABEL); print('-depsc',str)
    if (PREFS.PLOT_TCSPC == 1)
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
    str = sprintf('../output/%sTCSPC_chan.eps', PREFS.DESC_LABEL); print('-depsc',str)
    if (PREFS.PLOT_TCSPC == 1)
       figure;
    else
        clf;
    end;
end;

%
% =================================================================
%   MCS
% =================================================================

% MCS Plots with and without the time gate.

BIN_TIME_INC = PARAMS.BIN_TIME * 1E6; % this gives it to us in nsec.
% Apply time gates to get logical filters
DDFilter = (Dlc >= DDTimeGate(1)) & (Dlc <= DDTimeGate(2)); 
DAFilter = (Alc >= DATimeGate(1)) & (Alc <= DATimeGate(2));
AAFilter = (Alc >= AATimeGate(1)) & (Alc <= AATimeGate(2));
% ADfilter is meaningless

% these xXFilter arrays are over the list of photons.
% DDFilter is one when the photon is within the DDTimeGate range

% As a check we could do something like:
% (DDFilter | DAFilter | AAFilter) == ones()


% >> return maps that will be used for photon selection
DDFilterPos = find(DDFilter); % give me the indices of all photons that pass the DDFilter.
DAFilterPos = find(DAFilter);
AAFilterPos = find(AAFilter);
% AAFilterPos = filterPos(AAFilter);
% we have now assigned events

fprintf(1,'Number of DD counts: %i\n',length(DDFilterPos));
fprintf(1,'Number of DA counts: %i\n',length(DAFilterPos));
fprintf(1,'Number of AA counts: %i\n',length(AAFilterPos));

% >> global time
DDgt = Dgt(DDFilterPos); % pick 
DAgt = Agt(DAFilterPos);
AAgt = Agt(AAFilterPos);
% >> local time
DDlc = Dlc(DDFilterPos);
DAlc = Alc(DAFilterPos);
AAlc = Alc(AAFilterPos); %note that AA has a seriously right-shifted tcal.
% We have assigned the events to the proper channel.
% Can verify using (for example):
%   myHist = histc(DDlc,chan);
%   plot(chan,myHist);

% >> assign events in DD, DA,AA to an MCS bins.
maxTime = min([Dgt(end);Agt(end)]/1E6); % in milliseconds % divide by 60,000 to get minutes
fprintf(1,'Experiment time: %i min\n',floor(maxTime/60000));

% >> Calc the MCS data
[DDtmp, DDphotonMap, DDbinMapTmp] = calcMCS(DDgt, PARAMS.BIN_TIME, maxTime);
[DAtmp, DAphotonMap, DAbinMapTmp] = calcMCS(DAgt, PARAMS.BIN_TIME, maxTime);
[AAtmp, AAphotonMap, AAbinMapTmp]= calcMCS(AAgt, PARAMS.BIN_TIME, maxTime);

% A useful measure of noise is the number of ms with zero photons detected

% >> correct for possible differences in size of arrays XXtmp
smallest = min([length(DDtmp), length(DAtmp), length(AAtmp)]);
DD = DDtmp(1:smallest);
DA = DAtmp(1:smallest);
AA = AAtmp(1:smallest);
DDbinMap = DDbinMapTmp(1:smallest);
DAbinMap = DAbinMapTmp(1:smallest);
AAbinMap = AAbinMapTmp(1:smallest);
% This way, length(DD) = length(DA) = length(AA)

% NOTE, XXphotonMap is not corrected for "rogue" photons to an extra bin on
% the end.


meanCounts = [mean(DD), mean(DA), mean(AA)]; % gives average counts per bin
stdCounts = [std(DD), std(DA), std(AA)];

fprintf(1,'MCS statistics:');
display(meanCounts)
display(stdCounts)
diary off


if (PREFS.PLOT_BURST == 1) % plot individual burst traces?
    figure;
    start = 1000;
    len = 2000;
    str = sprintf('Counts / %1.1f ms',PARAMS.BIN_TIME);
    plotBurst1(PARAMS.BIN_TIME,DD,start,len,str); 
    title('DD');
    figure;
    plotBurst1(PARAMS.BIN_TIME,DA,start,len,str);
    title('DA');
    figure;
    plotBurst1(PARAMS.BIN_TIME,AA,start,len,str);
    title('AA');
end;

%
% =================================================================
%   PCH
% =================================================================

% ------------------ DD photon counting histogram ------------------

countBins = 0:max(DD); % Maximum y value
DDCounts = myHistc(DD,countBins)/maxTime;
DDCountsCum = cumsum(DDCounts);
semilogx(DDCounts,countBins,'b');
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
xlabel('Frequency [Events/sec]');
ylab = sprintf('Counts / %1.1f ms',PARAMS.BIN_TIME); ylabel(ylab); 
title('DD photon counting histogram');
grid on;
ax = axis;
hold on;
cMin = 1E-8;
cMax = max(DDCounts);
semilogx([cMin cMax], [PARAMS.NOISE_MIN PARAMS.NOISE_MIN],'r');
semilogx([cMin cMax], [PARAMS.NOISE_MAX PARAMS.NOISE_MAX],'g');
%xlim([cMin cMax]);
hold off;
axis(ax);
% print -depsc2 '../output/DDphotonCountingHistogram.eps'
str = sprintf('../output/%sDDphotonCountingHistogram.eps', PREFS.DESC_LABEL); print('-depsc',str);
if (PREFS.PLOT_PHOTON_COUNTING_HISTOGRAM == 1)
    figure;
else
    clf;
end;

% ------------------ DA photon counting histogram ------------------

countBins = 0:max(DA); % this produces max(DA)+1 bins.
DACounts = myHistc(DA,countBins)/maxTime;
semilogx(DACounts,countBins,'b');
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
xlabel('Frequency [Events/sec]');
ylab = sprintf('Counts / %1.1f ms',PARAMS.BIN_TIME); ylabel(ylab); 
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
str = sprintf('../output/%sDAphotonCountingHistogram.eps', PREFS.DESC_LABEL); print('-depsc',str);
if (PREFS.PLOT_PHOTON_COUNTING_HISTOGRAM == 1)
    figure;
else
    clf;
end;

% ------------------ DD + DA photon counting histogram ------------------

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
semilogx([cMin cMax], [PARAMS.NOISE_MIN PARAMS.NOISE_MIN],'r');
semilogx([cMin cMax], [PARAMS.NOISE_MAX PARAMS.NOISE_MAX],'g');
% xlim([cMin cMax]);
hold off;
%axis tight;
axis(ax);
% print -depsc2 '../output/DDDAphotonCountingHistogram.eps'
str = sprintf('../output/%sDDDAphotonCountingHistogram.eps', PREFS.DESC_LABEL); print('-depsc',str);
if (PREFS.PLOT_PHOTON_COUNTING_HISTOGRAM == 1)
    figure;
else
    clf;
end;

% ------------------ AA photon counting histogram ------------------

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
semilogx([cMin cMax], [PARAMS.PIE_MIN PARAMS.PIE_MIN],'r');
semilogx([cMin cMax], [PARAMS.PIE_MAX PARAMS.PIE_MAX],'g');
%xlim([cMin cMax]);
hold off;
axis(ax);
% print -depsc2 '../output/AAphotonCountingHistogram.eps'
str = sprintf('../output/%sAAphotonCountingHistogram.eps', PREFS.DESC_LABEL); print('-depsc',str);
if (PREFS.PLOT_PHOTON_COUNTING_HISTOGRAM == 1)
    figure;
else
    clf;
end;

%
% =================================================================
%   Main Analysis
% =================================================================

TEBins = linspace(PLOT_PREFS.TE_LIM(1),PLOT_PREFS.TE_LIM(2),PLOT_PREFS.NUM_BINS);
SBins = linspace(PLOT_PREFS.S_LIM(1),PLOT_PREFS.S_LIM(2),PLOT_PREFS.NUM_BINS);

% ------------------- Burst selection filters --------------------


% these are logical arrays that will be used to select bursts
redFilter = ((AA > PARAMS.PIE_MIN) & (AA < PARAMS.PIE_MAX)); % Acceptor is there
notRedFilter = (AA <= PARAMS.NOT_RED) & (DA <= PARAMS.NOT_RED); % identifies sample c/ no red counts
greenFilter = ((DD+DA > PARAMS.NOISE_MIN) & (DD+DA < PARAMS.NOISE_MAX)); % burst selection. Fluorescence could come out of donor or acceptor
notGreenFilter = (DD <= PARAMS.NOT_GREEN); % sample with no green counts

greenAndRedFilter = redFilter & greenFilter; % this gives the FRET-quenched donors
greenOrRedFilter = redFilter | greenFilter; % unused
greenAndNotRedFilter = greenFilter & notRedFilter; % this gives the donor-only population.
redAndNotGreenFilter = redFilter & notGreenFilter; % this gives acceptor-only population.

% >> return maps that will be used for photon selection
redFilterPos = find(redFilter);
greenFilterPos = find(greenFilter);
greenAndRedFilterPos = find(greenAndRedFilter);
greenOrRedFilterPos = find(greenOrRedFilter);
greenAndNotRedFilterPos = find(greenAndNotRedFilter);
% we have now assigned events

fprintf(1,'Total bursts: %i\n',length(greenFilterPos));
fprintf(1,'Total PIE selected bursts: %i\n',length(greenAndRedFilterPos));
fprintf(1,'Total Donor-only bursts: %i\n',length(greenAndNotRedFilterPos));



% ------------------- Burst Traces ----------------------------------

% >> without PIE filtering
%figure
plotBurst2sel(PARAMS.BIN_TIME,DD,DA,greenFilter,PLOT_PREFS);
pbaspect([4.0 1 1]);
title('Donor excitation'); %GHK
%print -depsc2 '../output/normalBurst.eps';
str = sprintf('../output/%sgreenExcitedBurst.eps', PREFS.DESC_LABEL); print('-depsc',str);
if (PREFS.PLOT_BURST_SEL == 1)
    figure;
else
    clf;
end;

% >> with PIE filtering
% Plot results (as subplots)
% top is burst plot of donor and acceptor with Donor excitation
% figure
subplot(2,1,1);
%plotBurst2sel(tstep,DD,DA,pieFilter,start,len);
plotBurst2sel(PARAMS.BIN_TIME,DD,DA,greenAndRedFilter,PLOT_PREFS);
title('Donor excitation');

% set(h,'XTick',[])
% set(gca,'XTickLabel',[]);

% bottom plot is plot of acceptor getting directly excited
subplot(2,1,2); 
plotBurst1sel(PARAMS.BIN_TIME,AA,redFilter,PLOT_PREFS);
title('Acceptor excitation');
% pbaspect([4.0 1 1]);

str = sprintf('../output/%spieBurst.eps', PREFS.DESC_LABEL); print('-depsc',str);
if (PREFS.PLOT_BURST_SEL == 1)
    figure;
else
    clf;
end;


% =================================================================
%   1D Histograms
% =================================================================


% >> 1D analysis with the PIE selection on.
F_D1 = PARAMS.GAMMA*DD + DA;
F_D2 = DD + DA; 
F_A = AA;

TE = DA./F_D1;
S = F_D2./(F_D2 + F_A); % S uses a gamma uncorrected DD.

% In other words:
% TE = DA./(gamma*DD + DA)
% S = (DD + DA)/(DD + DA + AA)

fTE = TE(greenFilterPos); % D-excited
fS = S(greenFilterPos);;
fTEPIE = TE(greenAndRedFilterPos); % D-excited and PIE selected
fSPIE = S(greenAndRedFilterPos);
fNTEPIE = TE(greenAndNotRedFilterPos); % DO-sample
fNSPIE = S(greenAndNotRedFilterPos);


TEHist = myHistc(fTE,TEBins); % noise-selected
SHist = myHistc(fS,SBins);

TEHistPIE = myHistc(fTEPIE,TEBins); %with PIE selection
SHistPIE = myHistc(fSPIE,SBins);

TEHistNotPIE = myHistc(fNTEPIE,TEBins); % DO-population
SHistNotPIE = myHistc(fNSPIE,SBins);

% assert sum(TEHistNotPIE) == sum(SHistNotPIE) == length(fNTEPIE) == T

% >> Plot of TE
figure(1),
subplot(2,1,1);
TEPlotPIE = bar(TEBins,TEHistPIE,1,'histc');
% >> plot the DO population
% hold on
% bar(TEBins,TEHistNotPIE,1,'histc');
% hole off
set(TEPlotPIE,'EdgeColor','k');
set(TEPlotPIE,'Facecolor','w');
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
%xlim([TEBins(1) TEBins(end)]);
xlim(PLOT_PREFS.TE_LIM);
xlabel('TE'); ylabel('Counts');
title('With PIE selection');

subplot(2,1,2);
TEPlot = bar(TEBins,TEHist,1,'histc');
set(TEPlot,'EdgeColor','k');
set(TEPlot,'Facecolor','w');
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
%xlim([TEBins(1) TEBins(end)]);
xlim(PLOT_PREFS.TE_LIM);
xlabel('TE'); ylabel('Counts');
title('Without PIE selection');
if PREFS.COMMENT
    legend(DESC);
end;
str = sprintf('../output/%s1D_TE.eps', PREFS.DESC_LABEL); print('-depsc',str);
if (PREFS.PLOT_1D_HIST == 1)
    figure;
else
    clf;
end;

% >> Plot of S
figure(2),
subplot(2,1,1);
SPlotPIE = bar(SBins,SHistPIE,1,'histc');
set(SPlotPIE,'EdgeColor','k');
set(SPlotPIE,'Facecolor','w');
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
%xlim([SBins(1) SBins(end)]);
xlim(PLOT_PREFS.S_LIM);
xlabel('S'); ylabel('Counts'); 
title('With PIE selection');

subplot(2,1,2);
SPlot = bar(SBins,SHist,1,'histc');
set(SPlot,'EdgeColor','k');
set(SPlot,'Facecolor','w');
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
%xlim([SBins(1) SBins(end)]);
xlim(PLOT_PREFS.S_LIM);
xlabel('S'); ylabel('Counts'); 
title('Without PIE selection');
if (PREFS.COMMENT == 1)
    legend(DESC);
end;
str = sprintf('../output/%s1D_S.eps', PREFS.DESC_LABEL); print('-depsc',str);
if PREFS.PLOT_1D_HIST
    figure;
else
    clf;
end;

%%
% =================================================================
%   2D Histograms
% =================================================================

% h = 2DHistPlot(@mainHistLog,TEBins,TEHist,SBins,SHist,fTE,fS,PREFS,...
%    '2D histogram (Log) w/o PIE selection','flipud(hot(256))');

% PLOT_PREFS.CMAP = flipud(hot(256));
% PLOT_PREFS.CMAP = flipud(pink(256));
% PLOT_PREFS.CMAP = flipud(bone(256));
PLOT_PREFS.CMAP = jet(256);
h2D = f2DHistPlot(@mainHist,TEBins,TEHist,SBins,SHist,fTE,fS,PREFS,PLOT_PREFS,'2D histogram','2DHist');

% PLOT_PREFS.CMAP = jet(256);
h2DL = f2DHistPlot(@mainHistLog,TEBins,TEHist,SBins,SHist,fTE,fS,PREFS,PLOT_PREFS,'2D histogram (Log)','2DHistLog');

% PLOT_PREFS.CMAP = jet(256);
h2DL = f2DHistPlot(@mainHist,TEBins,TEHistPIE,SBins,SHistPIE,fTEPIE,fSPIE,PREFS,PLOT_PREFS,'2D histogram with PIE selection','2DHistPIE');

% this is the DO-populaiton
h2DL = f2DHistPlot(@mainHist,TEBins,TEHistNotPIE,SBins,SHistNotPIE,fNTEPIE,fNSPIE,PREFS,PLOT_PREFS,'2D histogram zero peak selection','2DHistNotPIE');


%%
% =================================================================
%   Peak point selection
% =================================================================


% >> Peak point selection %GHK

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
xlim(PLOT_PREFS.TE_LIM);
xlabel('TE'); ylabel('Counts');
title(['TE of ',inFile.name]);
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

%%
% =================================================================
%   Fitting to 1 or more Gaussians
% =================================================================
%use matlab curvefitting toolbox to fit to gaussians, Feb. 09 2012
%
% Here, what we are doing is doing histogrms on center instead of on edges.
%
% PARAMETERS
figure;
NUM_GAUSSIANS = 1;     % 2 for two gaussian fit
REDUCED_FIT_RANGE = 1; % if true the don't fit to the full range of TE;s
%                      % 0 for full range, 1 for reduced fit range 
%
R_NUM_BINS = PLOT_PREFS.NUM_BINS - 1; % doing on centers reduces NUM_BINS by 1
binWidthD2 = (TEBins(2)-TEBins(1))/2;
TEMidBins = linspace(PLOT_PREFS.TE_LIM(1)+binWidthD2,PLOT_PREFS.TE_LIM(2)-binWidthD2,R_NUM_BINS);
%
if REDUCED_FIT_RANGE
%TEFitRange = [0.5 0.98]; % reduced fit range for 0.8 above
TEFitRange = [0.2 0.8]; % center weighted fitting
%TEFitRange = [0.02 0.55]; % reduced fit range for 0.3 below
else
   TEFitRange = PLOT_PREFS.TE_LIM;
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
title(['TE of ',inFile.name]);

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
% PARAMS.PIE_MIN, PARAMS.NOISE_MIN,meanTE
MTE = meanTE; % for texting mean TE on figure
ytext = max(TEHistPIE);
text(0.01,ytext,['Mean TE = ', num2str(MTE)],'FontSize',14,'BackgroundColor',[.7 .9 .7]);
% for print output folder
str = sprintf('../output/%sGaussianFittedMeanTE.eps', PREFS.DESC_LABEL); print('-depsc',str); %print on output folder
%
%


%%
% =================================================================
%   TCSPC of PIE-filtered data
% =================================================================

% greenAndRedFilterPos is an array of bins 

index = 1;
for i=1:length(greenAndRedFilterPos)
    m = greenAndRedFilterPos(i);
    mL = length(m);
    DDlcSel(index:(index+mL-1))=DDlc(m);
    index = index+mL;
end
