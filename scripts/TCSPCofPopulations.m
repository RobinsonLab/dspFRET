% =================================================================
%   TCSPC of PIE-filtered data
% =================================================================

% run this file after you have run analyzePIEgen2.m
% INPUTS:
%   things like: greenAndRedFilterPos, and various Maps.


% greenAndRedFilterPos is an array of bins 

IRF = importdata('../input/APD2_IRF.dat');

IRF1 = importdata('../input/IRF of APD2 (with 550).dat');
IRF2 = importdata('../input/IRF of APD2 (without 550).dat');
%%
% =================================================================
%   original DD, and sorted D and DA populations
% =================================================================

DDindex = 1;
DAindex = 1;
DOindex = 1;

DDlcSel = [];
DAlcSel = [];
DOlcSel = [];

for i=1:length(greenAndRedFilterPos)
    bQ = greenAndRedFilterPos(i); % bin of a FRET-quenched dample (DA sample)
    
    DDphotonsInBin = DDbinMap{bQ}; % get the photons that have mapped to the bin    
    DDl = length(DDphotonsInBin);
    DDlcSel(DDindex:(DDindex+DDl-1)) = DDlc(DDphotonsInBin);
    DDindex = DDindex + DDl;

    DAphotonsInBin = DAbinMap{bQ}; % get the photons that have mapped to the bin
    DAl = length(DAphotonsInBin);
    DAlcSel(DAindex:(DAindex+DAl-1)) = DAlc(DAphotonsInBin);
    DAindex = DAindex + DAl;
%    tot = DDl + DAl; % better be >= PARAMS.NOISE_MIN
end

for i=1:length(greenAndNotRedFilterPos)    
    bNQ = greenAndNotRedFilterPos(i); % bin of a non-FRET-quenched dample (DO sample)
    
    DOphotonsInBin = DDbinMap{bNQ}; % get the photons that have mapped to the bin
    DOl = length(DOphotonsInBin);
    DOlcSel(DOindex:(DOindex+DOl-1)) = DDlc(DOphotonsInBin);
    DOindex = DOindex + DOl;
end

t = chan*inFile.resolution; % in nsec
DDHist = histc(DDlcSel,chan); % do it on the center instead of edges (hist vs. histc)
DOHist = histc(DOlcSel,chan); % do it on the center instead of edges (hist vs. histc)
semilogy(t,donorHist,'k','linewidth',2);
hold on
semilogy(t,DDHist,'b','linewidth',2);
semilogy(t,DOHist,'g','linewidth',2);
hold off

% set(gcf,'linewidth',2);
set(gca,'FontName','Helvetica');
set(gca,'FontSize',18);
xlabel('Time [ns]');
ylabel('Counts');
grid on;
str = sprintf('../output/%sDD_DO_DA_TCSPC.eps', PREFS.DESC_LABEL); print('-depsc',str)
%%
% =================================================================
%   normalized on a log scale: a failure
% =================================================================

C = 10000;
DDHistN = C*DDHist/max(DDHist)+1;
DOHistN = C*DOHist/max(DOHist)+1;


% DDHistN = DDHist/max(DDHist);
% DOHistN = DOHist/max(DOHist);

numChan = length(chan);
%baseline = min(IRF(1:numChan-2));
baseline = mean(IRF(numChan/2:numChan-2));

IRFc = IRF - baseline;
min(IRFc);
IRFN = C*IRFc/max(IRF)+1;
% IRFN = IRFc/max(IRF);

semilogy(t,DDHistN,'b','linewidth',2);
hold on;
semilogy(t,DOHistN,'g','linewidth',2);
semilogy(t,IRFN(1:numChan),'r','linewidth',2);
hold off


%% 
% =================================================================
%   Linear scale, normalized
% =================================================================

C = 10000;
% DDHistN = C*DDHist/max(DDHist)+1;
% DOHistN = C*DOHist/max(DOHist)+1;

[DDHistM,DDHistI] = max(DDHist);
DDHistN = DDHist/DDHistM;

[DOHistM,DOHistI] = max(DOHist);
DOHistN = DOHist/DOHistM;

numChan = length(chan);
baseline = min(IRF(1:numChan-2));
% baseline = mean(IRF(numChan/2:numChan-2));

IRF = IRF2;

IRFc = IRF - baseline;
min(IRFc);
% IRFN = C*IRFc/max(IRF)+1;
[IRFM,IRFI] = max(IRF);
IRFN = IRFc/IRFM;

IRFI - DDHistI

% tcal= (IRFI - DDHistI)*inFile.resolution;
tcal= -21*inFile.resolution;

plot(t+tcal,DDHistN,'b','linewidth',2);
hold on;
plot(t+tcal,DOHistN,'g','linewidth',2);
plot(t,IRFN(1:numChan),'r--','linewidth',2);
hold off

xlim([1,6]);
ylim([0,1.1]);

set(gca,'FontName','Helvetica');
set(gca,'FontSize',18);
xlabel('Time [ns]');
ylabel('Normalized Counts');
grid on
str = sprintf('../output/%snormalizeTCSPC.eps', PREFS.DESC_LABEL); print('-depsc',str)



%%
% =================================================================
%   "Fits" log scale normalized
% =================================================================
tcal= -20*inFile.resolution;

semilogy(t+tcal,DDHistN,'b','linewidth',2);
hold on;
semilogy(t+tcal,DOHistN,'g','linewidth',2);
hold off;
xlim([1,6]);
ylim([0.003,1.1]);

% Donor-only
tau = 0.2;
y = exp(-t/tau);
fac = 4.5;

IRF1c = IRF1 - fac*min(IRF1(1:numChan-1));

decay = conv(IRF1c,y);

hold on;
semilogy(t,decay(1:numChan)/max(decay(1:numChan)),'k','linewidth',1.5);
hold off;

% Donor-acceptor
tau = 0.11;
y = exp(-t/tau);
fac = 12;
IRF1c = IRF1 - fac*min(IRF1(1:numChan-1));

decay = conv(IRF1c,y);

hold on;
semilogy(t,decay(1:numChan)/max(decay(1:numChan)),'k','linewidth',1.5);
hold off;

hold on;
semilogy(t,IRF1c(1:numChan)/max(IRF1c(1:numChan)),'r--','linewidth',2);
hold off;

xlim([1,6]);
ylim([0.01,1.1]);

set(gca,'FontName','Helvetica');
set(gca,'FontSize',18);
xlabel('Time [ns]');
ylabel('Normalized Counts');
grid on;
str = sprintf('../output/%snormalizeTCSPCFit.eps', PREFS.DESC_LABEL); print('-depsc',str)

%%
% =================================================================
%   IRF c/ and s/ the 550 clean up filter.
% =================================================================
semilogy(t,IRF1(1:numChan),'r','linewidth',2);
hold on;
semilogy(t,IRF2(1:numChan),'g','linewidth',2);
hold off;
xlim([1,6]);

%%

plot(t,IRF1(1:numChan)/max(IRF1(1:numChan)),'r','linewidth',2);
hold on;
plot(t,IRF2(1:numChan)/max(IRF2(1:numChan)),'g','linewidth',2);
hold off;
xlim([1,6]);

