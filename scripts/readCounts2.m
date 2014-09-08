function [Dlc Dgt Alc Agt] = readCounts2(PREFS, inFile);

% DESCRIPTION:
% This script reads a PicoHarp 300 T3 Mode data file (*.pt3)
% Works with file format version 2.0 only!

%
% REVISION HISTORY:
%
% Ver  Author(s)    Date   Description 
% ---  ----------   ----   --------------------------------------------
%  0    P. Kapusta  05/07  o original file name: read_pt3.m (PicoQuant GmbH)
%       M. Wahl
%  1    J. Robison  1/4/12 o Tested, speed up of at least 100x
%  2    J. Robison  6/8/12 o implemented data strcutures.
%                           o writing to a log file.


% NOTES:
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

% INPUTS:
% inFile- ptr, pointer to input file
%       - numRecords, number of records in input file to process
%       - resolution, timing resolution (in ns)
%       - syncPeriod, Sync period (in ns)

% PREFS- VERBOSE?

% OUTPUTS:
% 
% Dlc = local time (in channel number) of event in donor channel
% Dgt = global time (trueTime) of event in donor channel
% Alc = local time (in channel number) of event in acceptor channel
% Agt = global time (trueTime) of event in acceptor channel

% Why channel # instead of time?
% We express local time in channel # because we want to pick and filter in terms
% of channels

% Constants
WRAPAROUND=65536; % = 2^16 this is the maximum number of counts/sec that the PicoHarp can handle.
DIVISIONS = 10;

if PREFS.VERBOSE
   fprintf(1,'\n---------------------------------------------------------\n');
   fprintf(1,'>> readCounts2.m');
   fprintf(1,'\n---------------------------------------------------------\n');
   fprintf(1,'\nProcessing %d records\n',inFile.numRecords);
   fprintf(1,'Timeing Resolution: %5.6f ns\n', inFile.resolution);
   fprintf(1,'Sync Period = %5.4f ns\n',inFile.syncPeriod);
end

% Counters
% ofltime = 0;
% cnt_1=0; cnt_2=0; cnt_3=0; cnt_4=0; cnt_Ofl=0; cnt_M=0; cnt_Err=0;

% >> read in every single record from inFile into memory.
% T3Record is an array of size count of 32 bit integers

%     4.0265
%     4.0265
%     0.6177
%     0.5624
%     4.0265
%     0.5483
%     0.5458
%     4.0265
%     4.0265
%     4.0265


[T3Record count] = fread(inFile.ptr, inFile.numRecords, 'ubit32');     % all 32 bits:
% >> quality control
if count ~= inFile.numRecords
    fprintf(1,'\nFailed: could only read %d records of %d\n',count, inFile.numRecords);
    STOP;
end;

% now interpret T3Record
% everything is done on an array-wide basis.
% >> for each event, define the global time, the local time, and the channel.


%   +-------------------------------+  +-------------------------------+ 
%   |x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|  |x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|
%   +-------------------------------+  +-------------------------------+  
  
nsync = bitand(T3Record,65535);       % the lowest 16 bits:

% nsync gives which pulse just produced the event. 
% the start time of this pulse is nsync*syncPeriod.
% the problem is that nsync can only go up to 
  
%   +-------------------------------+  +-------------------------------+ 
%   | | | | | | | | | | | | | | | | |  |x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|x|
%   +-------------------------------+  +-------------------------------+  
  
chan = bitand(bitshift(T3Record,-28),15);   % the upper 4 bits:

%   +-------------------------------+  +-------------------------------+ 
%   |x|x|x|x| | | | | | | | | | | | |  | | | | | | | | | | | | | | | | |
%   +-------------------------------+  +-------------------------------+

%fprintf(outFile,'\n%7u %08x %6.0f  %2u   ',i,T3Record,nsync,chan);   

dChan = bitand(bitshift(T3Record,-16),4095);

% if chan = 1,2,3 or 4, then these  bits contain the dtime:
%   +-------------------------------+  +-------------------------------+ 
%   | | | | |x|x|x|x|x|x|x|x|x|x|x|x|  | | | | | | | | | | | | | | | | |
%   +-------------------------------+  +-------------------------------+    
    

% >> rollFilter, donFilter, and accFilter are masks 
rollFilter = (chan == 15); % this is a rollover
donFilter = (chan == 2); % count is from the donor channel
accFilter = (chan == 1); % count is from the acceptor channel
if all(rollFilter | donFilter | accFilter)
    fprintf(1,'\nAll %i records could be interpreted\n',count);
else
    fprintf(1,'\nOnly %i records could be interpreted\n',count);
    STOP;
end;

roll = transpose(zeros(1, inFile.numRecords));
roll(rollFilter) = WRAPAROUND;
cumRoll = cumsum(roll);
truenSync = cumRoll + nsync; % this is a count of how many times the laser has fired
% >> photon arrival event is the syncTime + dTime
% dTime = dChan*resolution;
trueTime = truenSync*inFile.syncPeriod + dChan*inFile.resolution; % in nsec  

Dlc = dChan(donFilter);
Dgt = trueTime(donFilter);
Alc = dChan(accFilter);
Agt = trueTime(accFilter);

if PREFS.VERBOSE
    fprintf(1,'All %i records have been processed\n',count);
end
 % MY NOTES
 % syncPeriod is ns/period
 % resolution is ns/chan
 % note that max(dtime)*resolution <= syncPeriod.
 % max(dtime) can not exceed the number of timing bins, which is 65536.
 % goal is to keep max(dtime) < 2k = 2048.
 % fprintf(outFile,' %10.0f %12.3f', truensync, truetime);
 

if PREFS.VERBOSE
    fprintf(1,'\nStatistics obtained from the data:\n');
%    fprintf(1,'\nLast True Sync = %-14.0f, Last t = %14.3f ns,',truenSync, trueTime);
    fprintf(1,'Ch1 (acceptor SPAD): %i counts, Ch2 (dondor SPAD): %i counts\n',length(Alc),length(Dlc));
%    fprintf(1,'%i overflows, %i markers, %i illegal events. Total: %i records read.\n',cnt_Ofl,cnt_M,cnt_Err,cnt_1+cnt_2+cnt_3+cnt_4+cnt_Ofl+cnt_M+cnt_Err);
%    fprintf(1,'\nRtCh1 Duty cycle: %-0.7f\n',cnt_1/truenSync);
%    fprintf(1,'\nRtCh2 Duty cycle: %-0.7f\n',cnt_2/truenSync);
%    fprintf(1,'\n');
end

end