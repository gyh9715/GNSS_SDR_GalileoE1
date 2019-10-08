function settings = initSettings()
%Functions initializes and saves settings. Settings can be edited inside of
%the function, updated from the command line or updated using a dedicated
%GUI - "setSettings".
%
%All settings are described inside function code.
%
%settings = initSettings()
%
%   Inputs: none
%
%   Outputs:
%       settings     - Receiver settings (a structure).

%--------------------------------------------------------------------------

%% Processing settings ====================================================
% Number of milliseconds to be processed used 32000 + any transients (see
% below - in Nav parameters) to ensure nav subframes are provided
settings.msToProcess        = 40000;        %[ms]      
% Number of channels to be used for signal processing
settings.numberOfChannels   = 6;

% Move the starting point of processing. Can be used to start the signal
% processing at any point in the data record (e.g. for long records). fseek
% function is used to move the file read point, therefore advance is byte
% based only. For Real sample files it skips the number of bytes as indicated
% here. For I/Q files it skips twice the number of bytes as indicated here
% to consider both I and Q samples
settings.skipTime              = 50; %unit s
settings.skipNumberOfBytes     = 10;

%% Raw signal file name and other parameter ===============================
% This is a "default" name of the data file (signal record) to be used in
% the post-processing mode  
    
% Data type used to store one sample
settings.dataType           = 'int8';
settings.fileName           = ...
   'D:\接收机资料\中频数据\L125_III1b_210s_L1.bin';
% File Types
%1 - 8 bit real samples S0,S1,S2,...
%2 - 8 bit I/Q samples I0,Q0,I1,Q1,I2,Q2,...
settings.fileType           = 2;

% Intermediate, sampling and code frequencies
settings.IF                 = 0e6;        %[Hz]
settings.carrFreqBasis      = 1575.42e6;  %[Hz]
settings.codeFreqBasis      = 1.023e6;    %[Hz]
settings.samplingFreq       = 20e6;       %[Hz]
settings.Datalength         = 1;          %[bytes]
settings.symbolRate         = int32 (250);%[bits/s]
 
% Define number of chips in a code period
settings.codeLength         = 4092;       %[chips]

%Calculate skiping samples
settings.skipNumberOfSamples= settings.samplingFreq * settings.skipTime;

%% Acquisition settings ===================================================
% Switch acquistion type
%0 - 8ms coherent integration time
%1 - 100ms coherent integration time
settings.acquisitionType    = 0;
% Skips acquisition in the script postProcessing.m if set to 1
settings.skipAcquisition    = 0;
% List of satellites to look for. Some satellites can be excluded to speed
% up acquisition
% settings.acqSatelliteList  
settings.acqSatelliteList   = 1:27;
% Band around IF to search for satellite signal. Depends on max Doppler
settings.acqSearchBand      = 14;           %[kHz]
% Threshold for the signal presence decision rule
settings.acqThreshold       = 2.5;

%% Tracking loops settings ================================================
% Skips tracking in the script postProcessing.m if set to 1
settings.SkipTracking    = 0;
% DLL code tracking loop parameters
settings.dllCorrelatorSpacing    =    2/12;     %[chips]
settings.dllVeryEarlyLateSpc     =    5/12;
settings.DLLa2                   =   1.414;
settings.DLLBL2                  = 0.74942;
% PLL carrier tracking loop parameters
settings.PLLa3                   =     1.1;
settings.PLLb3                   =     2.4;
settings.PLLBL3                  =    16.8;
% FLL carrier tracking loop parameters
settings.TrackingSettingFLLUsed  =  'true';
settings.FLLa2                   =   1.414;
settings.FLLBL2                  =       8;      %[Hz]

%% Navigation solution settings ===========================================

% Rate for calculating pseudorange and position
settings.navSolRate         = 5;            %[Hz]
settings.navSolPeriod       = 1000/settings.navSolRate; %ms

% Elevation mask to exclude signals from satellites at low elevation
settings.elevationMask      = 10;           %[degrees 0 - 90]
% Enable/dissable use of tropospheric correction
settings.useTropCorr        = 0;            
% 0 - Off
% 1 - On

% True position of the antenna in UTM system (if known). Otherwise enter
% all NaN's and mean position will be used as a reference .
settings.truePosition.E     = nan;
settings.truePosition.N     = nan;
settings.truePosition.U     = nan;

%% Plot settings ==========================================================
% Enable/disable plotting of the tracking results for each channel
settings.plotTracking       = 1;
% 0 - Off
% 1 - On

%% Constants ==============================================================
settings.c                  = 299792458;    % The speed of light, [m/s]
settings.startOffset        = 68.802;       %[ms] Initial sign. travel time
