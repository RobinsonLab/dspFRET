function figureHandle = f2DHistPlot(fHandle,TEBins,TEHist,SBins,SHist,fTE,fS,PREFS,PLOT_PREFS,gTitle,fName)


% INPUT: 
% fHandle, wither mainHist or mainHistLog

% Notes:
% for 2D histogram analysis, we select only based on countFilter

figureHandle = figure() % brings figure 4 to the front.

% >> TEHistPlot
subplot(2,2,1); 
bar(TEBins,TEHist,1,'histc'); 
TEHistPlot = gca; 
axis([PLOT_PREFS.TE_LIM 0 max(TEHist)*1.01]); %axis('off');
%set(TEHistPlot,'xtick',[]); % don't write the x axis
set(gca,'FontName','Helvetica');
set(gca,'FontSize',14);
ylabel('Counts');
title(gTitle);
% grid on;

% >> Main plot
subplot(2,2,3); 
%mainHist(fTE,fS,TEBins,SBins,1,'k.',jet(256));
%mainHist(fTE,fS,TEBins,SBins,1,'k.',flipud(bone(256)));
fHandle(fTE,fS,TEBins,SBins,1,'k.',PLOT_PREFS.CMAP);  % Hot color
h1 = gca; % axis([xlim ylim]);
xlabel('E'); ylabel('S');

% >> SHistPlot
subplot(2,2,4);     
barh(SBins,SHist,1,'histc');
SHistPlot = gca; 
axis([0 max(SHist)*1.01 PLOT_PREFS.S_LIM]); %axis('off');
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

str = sprintf('../output/%s%s.eps', PREFS.DESC_LABEL,fName); print('-depsc2',str);
% str = sprintf('../output/%s%s.jpg', PREFS.DESC_LABEL,fName); print('-djpeg90',str);
% open in AI -> RBG. Then export as 300 dpi png

if PREFS.COMMENT
    bb = legend(PREFS.DESC_LABEL);
    set (bb, 'Position',[.70 .82 .20 .05]);
end;
%str = sprintf('../output/%s2DHistThresh.eps', DESC_LABEL); print('-depsc',str);

