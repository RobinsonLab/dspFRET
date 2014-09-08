% FILENAME: processpt3.m
%
% DESCRIPTION: This is the master file for reading in and processing a .pt3
% file
%
% DEPENDENT FILES: 
%   readHeader
%   readCounts2
%
%
% INPUT FILES
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
clc;

% CONSTANTS

% BEGIN


% USER VARIABLES

GUI_PREF = 1 % 1 = true; 0 = false
% Set filename and pathname using the GUI

if (GUI_PREF == 1)
    [filename, pathname]=uigetfile('*.pt3', 'T3 Mode data:');
else
% Set filename and pathname manually
    filename = 'Annealing buffer(5%glycerol).pt3';
    pathname = '../input/';
end;
inFile=fopen([pathname filename]);
logFileName = [pathname filename(1:length(filename)-4) '.log'];
logFile = fopen(logFileName,'W');
fprintf(logFile,'Header from file: %s\n', filename);
fprintf(logFile,'\n');
[numRecords resolution syncPeriod] = readHeader(inFile, logFile);
[Ddt Dtt Adt Att] = readCounts2(inFile, logFile, numRecords, resolution, syncPeriod);
fclose(inFile);
fclose(logFile);

% Analyze the Background

% % if you want in terms of time [nsec]
% chan = 0:resolution:syncPeriod;
% donorHist = histc(Ddt*resolution,chan);
% acceptorHist = histc(Adt*resolution,chan);

FWHM = 7; 
% TCSPC plots in terms of channel #
% allows us to the channel # for time gating
chan = 0:round(syncPeriod/resolution); % gives integers
t = chan*resolution; % in nsec
donorHist = histc(Ddt,chan); % note: Ddt is in channels, not time.
acceptorHist = histc(Adt,chan);
% USER VARIABLES
myOffset = 700;
[maxRng,maxRngPos] = max(acceptorHist(myOffset:1200)); 
PIEPeak = myOffset + maxRngPos;
AATimeGate = [PIEPeak + FWHM, length(chan)-10];
[maxRng,maxRngPos] = max(donorHist(1:600));
DDTimeGate = [maxRngPos + FWHM, PIEPeak - 3*FWHM];
[maxRng,maxRngPos] = max(acceptorHist(1:600));
DATimeGate = [maxRngPos + FWHM, PIEPeak - 3*FWHM];

spadShift = DDTimeGate(1) - DATimeGate(1); % just for fun

% green excitation
subplot(2,1,1);
rectangle('Position',[DDTimeGate(1),0,DDTimeGate(2)-DDTimeGate(1),max(donorHist)*1.1],'FaceColor',[0.9 0.9 0.9],'edgecolor','none');
hold on;
plot(chan,donorHist,'b');
% semilogy(timeChan,donorHist,'b');
hold off;
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
ylabel('Counts'); xlabel('Channel'); 
title('Donor');
grid on;

subplot(2,1,2);
rectangle('Position',[DATimeGate(1),0,DATimeGate(2)-DATimeGate(1),max(acceptorHist)*1.1],'FaceColor',[0.9 0.9 0.9],'edgecolor','none');
hold on;
rectangle('Position',[AATimeGate(1),0,AATimeGate(2)-AATimeGate(1),max(acceptorHist)*1.1],'FaceColor',[0.9 0.9 0.9],'edgecolor','none');
plot(chan,acceptorHist,'r');
hold off;
% semilogy(timeChan,acceptorHist,'r');
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
ylabel('Counts'); xlabel('Channel'); 
title('Acceptor');
grid on;
print -depsc2 '../output/TCSPCchan.eps'

% OUTPUTS

% DDTimeGate, DATimeGate, AATimeGate
save '../input/timeGateFilters.mat' DDTimeGate DATimeGate AATimeGate;
whos('-file', '../input/timeGateFilters.mat') %check that they were written
% % plot in time works
% figure
% plot(t,donorHist,'b');
% set(gca,'FontName','Helvetica');
% set(gca,'FontSize',14);
% ylabel('Counts'); xlabel('Time [nsec]'); 
% title('Donor');

%% Analyze the data.

if (GUI_PREF == 1)
    [filename, pathname]=uigetfile('*.pt3', 'T3 Mode data:');
else
% Set filename and pathname manually
    filename = 'Annealing buffer(5%glycerol).pt3';
    pathname = '../input/';
end;
inFile=fopen([pathname filename]);
logFileName = [pathname filename(1:length(filename)-4) '.log'];
logFile = fopen(logFileName,'W');
fprintf(logFile,'Header from file: %s\n', filename);
fprintf(logFile,'\n');
[numRecords resolution syncPeriod] = readHeader(inFile, logFile);
[Ddt Dtt Adt Att] = readCounts2(inFile, logFile, numRecords, resolution, syncPeriod);
fclose(inFile);
fclose(logFile);



%% MCS Plots with and without the time gate.

BIN_TIME = 10.0; % in millisecond
BIN_TIME_INC = BIN_TIME * 1E6; % this gives it to us in nsec.
% Apply time gates
DDFilter = (Ddt >= DDTimeGate(1)) & (Ddt <= DDTimeGate(2)); 
DAFilter = (Adt >= DATimeGate(1)) & (Adt <= DATimeGate(2));
AAFilter = (Adt >= AATimeGate(1)) & (Adt <= AATimeGate(2));

maxTime = min([Dtt(end);Att(end)]/1E6); % in milliseconds % divide by 60,000 to get minutes

DDcounts = getMCS(Dtt,DDFilter, BIN_TIME,maxTime);
DAcounts = getMCS(Att,DAFilter, BIN_TIME,maxTime);
AAcounts = getMCS(Att,AAFilter, BIN_TIME,maxTime);
meanCounts = [mean(DDcounts), mean(DAcounts), mean(AAcounts)] % gives average counts per bin
stdCounts = [std(DDcounts), std(DAcounts), std(AAcounts)]

str = sprintf('Counts / %d ms',BIN_TIME);
plotBurst1(BIN_TIME,DDcounts,1,5000,str);
title('DD');
figure;
plotBurst1(BIN_TIME,DAcounts,1,5000,str);
title('DA');
figure;
plotBurst1(BIN_TIME,AAcounts,1,5000,str);
title('DD');

%% Photon counting histogram
countBins = 0:max(DDcounts);
DDhist = myHistc(DDcounts,countBins)/maxTime;
figure
semilogx(DDhist,countBins,'b');
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
xlabel('Frequency (Events/sec)');
ylabel(str);
title('DD photon counting histogram');
grid on;
hold on;
mn = mean(DDcounts);
sd1 = mean(DDcounts) + 1*std(DDcounts);
sd2 = mean(DDcounts) + 2*std(DDcounts);
sd3 = mean(DDcounts) + 3*std(DDcounts);
sd4 = mean(DDcounts) + 4*std(DDcounts);
sd5 = mean(DDcounts) + 5*std(DDcounts);
cMin = 1E-6;
cMax = max(DDhist);
semilogx([cMin cMax], [mn mn],'-k');
semilogx([cMin cMax], [sd1 sd1],'-.k');
semilogx([cMin cMax], [sd2 sd2],'--k');
semilogx([cMin cMax], [sd3 sd3],':k');
semilogx([cMin cMax], [sd4 sd4],':g');
semilogx([cMin cMax], [sd5 sd5],':r');

% xlim([cMin cMax]);
hold off;
axis tight;

% THE END RESULT

% We know
% (1) time gates for DD, DA, and AA
% (2) MCS data for DD, DA, and AA
% (3) photon counting histogram of DD, DA, AA


% Question is:

% When we time gate, we lower our burst intensity but increase the S/N ratio. But by cutting down our intensity, this introduces noise in
% the TE histogram.

% The alternative is not to time gate. throwing away bursts?



%% Now we process the DA data.