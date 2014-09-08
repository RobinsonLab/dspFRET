function [Records Resolution syncPeriod] = readHeader(inFile,logFile)
% This function reads the header of a Symphotime TTTR binary file.

Ident = char(fread(inFile, 16, 'char'));
fprintf(logFile,'      Identifier: %s\n', Ident);

FormatVersion = deblank(char(fread(inFile, 6, 'char')'));
fprintf(logFile,'  Format Version: %s\n', FormatVersion);

if not(strcmp(FormatVersion,'2.0'))
   fprintf(logFile,'\n\n      Warning: This program is for version 2.0 only. Aborted.');
   return;
end;

CreatorName = char(fread(inFile, 18, 'char'));
fprintf(logFile,'    Creator Name: %s\n', CreatorName);

CreatorVersion = char(fread(inFile, 12, 'char'));
fprintf(logFile,' Creator Version: %s\n', CreatorVersion);

FileTime = char(fread(inFile, 18, 'char'));
fprintf(logFile,'       File Time: %s\n', FileTime);

CRLF = char(fread(inFile, 2, 'char'));

CommentField = char(fread(inFile, 256, 'char'));
fprintf(logFile,'         Comment: %s\n', CommentField);


%
% The following is binary file header information
%


Curves = fread(inFile, 1, 'int32');
fprintf(logFile,'Number of Curves: %d\n', Curves);

BitsPerRecord = fread(inFile, 1, 'int32');
fprintf(logFile,'   Bits / Record: %d\n', BitsPerRecord);

RoutingChannels = fread(inFile, 1, 'int32');
fprintf(logFile,'Routing Channels: %d\n', RoutingChannels);

NumberOfBoards = fread(inFile, 1, 'int32');
fprintf(logFile,'Number of Boards: %d\n', NumberOfBoards);

ActiveCurve = fread(inFile, 1, 'int32');
fprintf(logFile,'    Active Curve: %d\n', ActiveCurve);

MeasurementMode = fread(inFile, 1, 'int32');
fprintf(logFile,'Measurement Mode: %d\n', MeasurementMode);

SubMode = fread(inFile, 1, 'int32');
fprintf(logFile,'        Sub-Mode: %d\n', SubMode);

RangeNo = fread(inFile, 1, 'int32');
fprintf(logFile,'       Range No.: %d\n', RangeNo);

Offset = fread(inFile, 1, 'int32');
fprintf(logFile,'          Offset: %d ns \n', Offset);

AcquisitionTime = fread(inFile, 1, 'int32');
fprintf(logFile,'Acquisition Time: %d ms \n', AcquisitionTime);

StopAt = fread(inFile, 1, 'int32');
fprintf(logFile,'         Stop At: %d counts \n', StopAt);

StopOnOvfl = fread(inFile, 1, 'int32');
fprintf(logFile,'Stop on Overflow: %d\n', StopOnOvfl);

Restart = fread(inFile, 1, 'int32');
fprintf(logFile,'         Restart: %d\n', Restart);

DispLinLog = fread(inFile, 1, 'int32');
fprintf(logFile,' Display Lin/Log: %d\n', DispLinLog);

DispTimeFrom = fread(inFile, 1, 'int32');
fprintf(logFile,' Display Time Axis From: %d ns \n', DispTimeFrom);

DispTimeTo = fread(inFile, 1, 'int32');
fprintf(logFile,'   Display Time Axis To: %d ns \n', DispTimeTo);

DispCountFrom = fread(inFile, 1, 'int32');
fprintf(logFile,'Display Count Axis From: %d\n', DispCountFrom); 

DispCountTo = fread(inFile, 1, 'int32');
fprintf(logFile,'  Display Count Axis To: %d\n', DispCountTo);

for i = 1:8
DispCurveMapTo(i) = fread(inFile, 1, 'int32');
DispCurveShow(i) = fread(inFile, 1, 'int32');
end;

for i = 1:3
ParamStart(i) = fread(inFile, 1, 'float');
ParamStep(i) = fread(inFile, 1, 'float');
ParamEnd(i) = fread(inFile, 1, 'float');
end;

RepeatMode = fread(inFile, 1, 'int32');
fprintf(logFile,'        Repeat Mode: %d\n', RepeatMode);

RepeatsPerCurve = fread(inFile, 1, 'int32');
fprintf(logFile,'     Repeat / Curve: %d\n', RepeatsPerCurve);

RepeatTime = fread(inFile, 1, 'int32');
fprintf(logFile,'        Repeat Time: %d\n', RepeatTime);

RepeatWait = fread(inFile, 1, 'int32');
fprintf(logFile,'   Repeat Wait Time: %d\n', RepeatWait);

ScriptName = char(fread(inFile, 20, 'char'));
fprintf(logFile,'        Script Name: %s\n', ScriptName);


%
% The next is a board specific header
%


HardwareIdent = char(fread(inFile, 16, 'char'));
fprintf(logFile,'Hardware Identifier: %s\n', HardwareIdent);

HardwareVersion = char(fread(inFile, 8, 'char'));
fprintf(logFile,'   Hardware Version: %s\n', HardwareVersion);

HardwareSerial = fread(inFile, 1, 'int32');
fprintf(logFile,'   HW Serial Number: %d\n', HardwareSerial);

SyncDivider = fread(inFile, 1, 'int32');
fprintf(logFile,'       Sync Divider: %d\n', SyncDivider);

CFDZeroCross0 = fread(inFile, 1, 'int32');
fprintf(logFile,'CFD ZeroCross (Ch0): %4i mV\n', CFDZeroCross0);

CFDLevel0 = fread(inFile, 1, 'int32');
fprintf(logFile,'CFD Discr     (Ch0): %4i mV\n', CFDLevel0);

CFDZeroCross1 = fread(inFile, 1, 'int32');
fprintf(logFile,'CFD ZeroCross (Ch1): %4i mV\n', CFDZeroCross1);

CFDLevel1 = fread(inFile, 1, 'int32');
fprintf(logFile,'CFD Discr     (Ch1): %4i mV\n', CFDLevel1);

Resolution = fread(inFile, 1, 'float');
fprintf(logFile,'         Resolution: %5.6f ns\n', Resolution);

% below is new in format version 2.0

RouterModelCode      = fread(inFile, 1, 'int32');
RouterEnabled        = fread(inFile, 1, 'int32');

% Router Ch1
RtChan1_InputType    = fread(inFile, 1, 'int32');
RtChan1_InputLevel   = fread(inFile, 1, 'int32');
RtChan1_InputEdge    = fread(inFile, 1, 'int32');
RtChan1_CFDPresent   = fread(inFile, 1, 'int32');
RtChan1_CFDLevel     = fread(inFile, 1, 'int32');
RtChan1_CFDZeroCross = fread(inFile, 1, 'int32');
% Router Ch2
RtChan2_InputType    = fread(inFile, 1, 'int32');
RtChan2_InputLevel   = fread(inFile, 1, 'int32');
RtChan2_InputEdge    = fread(inFile, 1, 'int32');
RtChan2_CFDPresent   = fread(inFile, 1, 'int32');
RtChan2_CFDLevel     = fread(inFile, 1, 'int32');
RtChan2_CFDZeroCross = fread(inFile, 1, 'int32');
% Router Ch3
RtChan3_InputType    = fread(inFile, 1, 'int32');
RtChan3_InputLevel   = fread(inFile, 1, 'int32');
RtChan3_InputEdge    = fread(inFile, 1, 'int32');
RtChan3_CFDPresent   = fread(inFile, 1, 'int32');
RtChan3_CFDLevel     = fread(inFile, 1, 'int32');
RtChan3_CFDZeroCross = fread(inFile, 1, 'int32');
% Router Ch4
RtChan4_InputType    = fread(inFile, 1, 'int32');
RtChan4_InputLevel   = fread(inFile, 1, 'int32');
RtChan4_InputEdge    = fread(inFile, 1, 'int32');
RtChan4_CFDPresent   = fread(inFile, 1, 'int32');
RtChan4_CFDLevel     = fread(inFile, 1, 'int32');
RtChan4_CFDZeroCross = fread(inFile, 1, 'int32');

% Router settings are meaningful only for an existing router:

if RouterModelCode>0

    fprintf(logFile,'-------------------------------------\n'); 
    fprintf(logFile,'   Router Model Code: %d \n', RouterModelCode);
    fprintf(logFile,'      Router Enabled: %d \n', RouterEnabled);
    fprintf(logFile,'-------------------------------------\n'); 
    
    % Router Ch1 
    fprintf(logFile,'RtChan1 InputType   : %d \n', RtChan1_InputType);
    fprintf(logFile,'RtChan1 InputLevel  : %4i mV\n', RtChan1_InputLevel);
    fprintf(logFile,'RtChan1 InputEdge   : %d \n', RtChan1_InputEdge);
    fprintf(logFile,'RtChan1 CFDPresent  : %d \n', RtChan1_CFDPresent);
    fprintf(logFile,'RtChan1 CFDLevel    : %4i mV\n', RtChan1_CFDLevel);
    fprintf(logFile,'RtChan1 CFDZeroCross: %4i mV\n', RtChan1_CFDZeroCross);
    fprintf(logFile,'-------------------------------------\n'); 

    % Router Ch2
    fprintf(logFile,'RtChan2 InputType   : %d \n', RtChan2_InputType);
    fprintf(logFile,'RtChan2 InputLevel  : %4i mV\n', RtChan2_InputLevel);
    fprintf(logFile,'RtChan2 InputEdge   : %d \n', RtChan2_InputEdge);
    fprintf(logFile,'RtChan2 CFDPresent  : %d \n', RtChan2_CFDPresent);
    fprintf(logFile,'RtChan2 CFDLevel    : %4i mV\n', RtChan2_CFDLevel);
    fprintf(logFile,'RtChan2 CFDZeroCross: %4i mV\n', RtChan2_CFDZeroCross);
    fprintf(logFile,'-------------------------------------\n'); 

    % Router Ch3
    fprintf(logFile,'RtChan3 InputType   : %d \n', RtChan3_InputType);
    fprintf(logFile,'RtChan3 InputLevel  : %4i mV\n', RtChan3_InputLevel);
    fprintf(logFile,'RtChan3 InputEdge   : %d \n', RtChan3_InputEdge);
    fprintf(logFile,'RtChan3 CFDPresent  : %d \n', RtChan3_CFDPresent);
    fprintf(logFile,'RtChan3 CFDLevel    : %4i mV\n', RtChan3_CFDLevel);
    fprintf(logFile,'RtChan3 CFDZeroCross: %4i mV\n', RtChan3_CFDZeroCross);
    fprintf(logFile,'-------------------------------------\n'); 

    % Router Ch4
    fprintf(logFile,'RtChan4 InputType   : %d \n', RtChan4_InputType);
    fprintf(logFile,'RtChan4 InputLevel  : %4i mV\n', RtChan4_InputLevel);
    fprintf(logFile,'RtChan4 InputEdge   : %d \n', RtChan4_InputEdge);
    fprintf(logFile,'RtChan4 CFDPresent  : %d \n', RtChan4_CFDPresent);
    fprintf(logFile,'RtChan4 CFDLevel    : %4i mV\n', RtChan4_CFDLevel);
    fprintf(logFile,'RtChan4 CFDZeroCross: %4i mV\n', RtChan4_CFDZeroCross);
    fprintf(logFile,'-------------------------------------\n'); 
 
end;
 
%
% The next is a T3 mode specific header
%

ExtDevices = fread(inFile, 1, 'int32');
fprintf(logFile,'   External Devices: %d\n', ExtDevices);

Reserved1 = fread(inFile, 1, 'int32');
fprintf(logFile,'          Reserved1: %d\n', Reserved1);

Reserved2 = fread(inFile, 1, 'int32');
fprintf(logFile,'          Reserved2: %d\n', Reserved2);

CntRate0 = fread(inFile, 1, 'int32');
fprintf(logFile,'   Count Rate (Ch0): %d Hz\n', CntRate0);

CntRate1 = fread(inFile, 1, 'int32');
fprintf(logFile,'   Count Rate (Ch1): %d Hz\n', CntRate1);

StopAfter = fread(inFile, 1, 'int32');
fprintf(logFile,'         Stop After: %d ms \n', StopAfter);

StopReason = fread(inFile, 1, 'int32');
fprintf(logFile,'        Stop Reason: %d\n', StopReason);

Records = fread(inFile, 1, 'uint32');
fprintf(logFile,'  Number Of Records: %d\n', Records);

ImgHdrSize = fread(inFile, 1, 'int32');
fprintf(logFile,'Imaging Header Size: %d bytes\n', ImgHdrSize);

%Special header for imaging 
ImgHdr = fread(inFile, ImgHdrSize, 'int32');

syncPeriod = 1E9/CntRate0;   % in nanoseconds
fprintf(logFile,'Sync Rate = %d / second\n',CntRate0);
fprintf(logFile,'Sync Period = %5.4f ns\n',syncPeriod);


end


