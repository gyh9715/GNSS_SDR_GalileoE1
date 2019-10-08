function [trackResults, channel]= tracking(fid, channel, settings)
% Performs code and carrier tracking for all channels.
%[trackResults, channel] = tracking(fid, channel, settings)
%
%   Inputs:
%       fid             - file identifier of the signal record for I
%       channel         - PRN, carrier frequencies and code phases of all
%                       satellites to be tracked (prepared by preRum.m from
%                       acquisition results).
%       settings        - receiver settings.
%   Outputs:
%       trackResults    - tracking results (structure array). Contains
%                       in-phase prompt outputs and absolute spreading
%                       code's starting positions, together with other
%                       observation data from the tracking loops. All are
%                       saved every millisecond.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%
% Copyright (C) Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
% Based on code by DMAkos Oct-1999
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: tracking.m,v 1.14.2.31 2006/08/14 11:38:22 dpl Exp $

%% Initialize
% Total tracking counts
trackTimes = settings.msToProcess / 4;     % For Galileo one raning code period is 4 ms

% Channel status
trackResults.status         = '-';      % No tracked signal, or lost lock

% The absolute sample in the record of the C/A code start:
trackResults.absoluteSample = zeros(1, trackTimes);

% Freq of the C/A code:
trackResults.codeFreq       = inf(1, trackTimes);
trackResults.remCodePhase   = inf(1, trackTimes);%record codephase zsh

% Frequency of the tracked carrier wave:
trackResults.carrFreq       = inf(1, trackTimes);
trackResults.remCarrPhase   = inf(1, trackTimes);%record carrier phase zsh

% Outputs from the correlators (In-phase):
trackResults.I_P            = zeros(1, trackTimes);
trackResults.I_E            = zeros(1, trackTimes);
trackResults.I_L            = zeros(1, trackTimes);

% Outputs from the correlators (Quadrature-phase):
trackResults.Q_E            = zeros(1, trackTimes);
trackResults.Q_P            = zeros(1, trackTimes);
trackResults.Q_L            = zeros(1, trackTimes);

% Loop discriminators
trackResults.dllDiscr       = inf(1, trackTimes);
trackResults.dllDiscrFilt   = inf(1, trackTimes);
trackResults.pllDiscr       = inf(1, trackTimes);
trackResults.pllDiscrFilt   = inf(1, trackTimes);

% Initialize filter parameters
for channelNr = 1:settings.numberOfChannels
    
        %FLL code tracking loop parameters
        Tracking.trackResults(channelNr).trackingweOldNco_1=0;
        Tracking.trackResults(channelNr).trackingweOldNco_2=0;
        Tracking.trackResults(channelNr).trackingweOldError_1=0;
        Tracking.trackResults(channelNr).trackingweOldError_2=0;

        Tracking.trackResults(channelNr).I_P_1=0;
        Tracking.trackResults(channelNr).Q_P_1=0;
        
        % PLL carrier/Costas loop parameters
        Tracking.trackResults(channelNr).trackingcarrOldNco_1=0;
        Tracking.trackResults(channelNr).trackingcarrOldNco_2=0;
        Tracking.trackResults(channelNr).trackingcarrOldError_1=0;
        Tracking.trackResults(channelNr).trackingcarrOldError_2=0;
        
        %DLL code tracking loop parameters
        Tracking.trackResults(channelNr).trackingcodeOldNco_1=0;
        Tracking.trackResults(channelNr).trackingcodeOldError_1=0;
   
end


%--- Copy initial settings for all channels -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Initialize tracking variables ==========================================

%--- DLL variables --------------------------------------------------------
% Define early-late offset (in chips)
earlyLateSpc     = settings.dllCorrelatorSpacing;
veryEarlyLateSpc = settings.dllVeryEarlyLateSpc;


%--- PLL variables --------------------------------------------------------
% Coherent integration period
PDIcarr = 0.004; % unit:second 

% hwb = waitbar(0,'Tracking...','Visible','off');
hwb = waitbar(0,'Tracking...');

if (settings.fileType==1)
    dataAdaptCoeff=1;
else
    dataAdaptCoeff=2;
end

%% Start processing channels ==============================================
for channelNr = 1:settings.numberOfChannels

    % Only process if PRN is non zero (acquisition was successful)
    if (channel(channelNr).PRN ~= 0)
        % Save additional information - each channel's tracked PRN
        trackResults(channelNr).PRN     = channel(channelNr).PRN;

        % Move the starting point of processing. Can be used to start the
        % signal processing at any point in the data record (e.g. for long
        % records). In addition skip through that data file to start at the
        % appropriate sample (corresponding to code phase). Assumes sample
        % type is schar (or 1 byte per sample)
        fseek(fid, ...
              dataAdaptCoeff*settings.Datalength*(settings.skipNumberOfSamples + channel(channelNr).codePhase-1), ...
              'bof');



        % Get a vector with the PN code sampled 1x/chip
        dataChannelBocCode   = generatePNcode(channel(channelNr).PRN,1,1);
        pilotChannelCode     = generatePNcode(channel(channelNr).PRN,2,1);
        dataChannelBocCode   = [dataChannelBocCode(4092*2*6-11:4092*2*6),...
                                dataChannelBocCode,...
                                dataChannelBocCode(1:12)];
        pilotChannelCode  = [pilotChannelCode(4092*2*6-11:4092*2*6) ...
                                pilotChannelCode...
                                pilotChannelCode(1:12)];

        %--- Perform various initializations ------------------------------
        
        %initialize carrier setting
        Tracking.trackResults(channelNr).carrFreq = settings.IF;
        Tracking.trackResults(channelNr).carrierDoppler = channel(channelNr).acquiredFreq;
        Tracking.trackResults(channelNr).remCarrPhase = 0;
        
        %initialize code Freq setting
        Tracking.trackResults(channelNr).codeDoppler = Tracking.trackResults(channelNr).carrierDoppler * ...
            12 * settings.codeFreqBasis / settings.carrFreqBasis;  %注意这里乘以12  
        Tracking.trackResults(channelNr).remCodePhase  = 0;
        Tracking.trackResults(channelNr).I_P_1=0;
        Tracking.trackResults(channelNr).Q_P_1=0;
        
        %% === Process the number of specified code periods =================
        for loopCnt =  1:trackTimes
        %FLL carrier tracking loop parameters
        weOldNco_1    = Tracking.trackResults(channelNr).trackingweOldNco_1;
        weOldNco_2    = Tracking.trackResults(channelNr).trackingweOldNco_2;
        weOldError_1  = Tracking.trackResults(channelNr).trackingweOldError_1;
        weOldError_2  = Tracking.trackResults(channelNr).trackingweOldError_2;
        
        % PLL carrier/Costas loop parameters
        carrOldNco_1   = Tracking.trackResults(channelNr).trackingcarrOldNco_1;
        carrOldNco_2   = Tracking.trackResults(channelNr).trackingcarrOldNco_2;
        carrOldError_1 = Tracking.trackResults(channelNr).trackingcarrOldError_1;
        carrOldError_2 = Tracking.trackResults(channelNr).trackingcarrOldError_2;
        
        %DLL code tracking loop parameters
        codeOldNco_1   = Tracking.trackResults(channelNr).trackingcodeOldNco_1;
        codeOldError_1 = Tracking.trackResults(channelNr).trackingcodeOldError_1;
        
        I_P_1=Tracking.trackResults(channelNr).I_P_1;
        Q_P_1=Tracking.trackResults(channelNr).Q_P_1;

            %% GUI update -------------------------------------------------------------
            % The GUI is updated every 50ms. This way Matlab GUI is still
            % responsive enough. At the same time Matlab is not occupied
            % all the time with GUI task.
            if (rem(loopCnt, 100) == 0)
                try
                    waitbar(loopCnt/trackTimes, ...
                            hwb, ...
                            ['Tracking: Ch ', int2str(channelNr), ...
                            ' of ', int2str(settings.numberOfChannels), ...
                            '; PRN#', int2str(channel(channelNr).PRN), ...
                            '; Completed ',int2str(loopCnt*4), ...
                            ' of ', int2str(trackTimes*4), ' msec']);                       
                catch
                    % The progress bar was closed. It is used as a signal
                    % to stop, "cancel" processing. Exit.
                    disp('Progress bar closed, exiting...');
                    return
                end
            end

            %% Read next block of data ------------------------------------------------
            % Find the size of a "block" or code period in whole samples

            % Update the phasestep based on code freq (variable) and
            % sampling frequency (fixed)
            codeFreq = Tracking.trackResults(channelNr).codeDoppler + 12*settings.codeFreqBasis;%子载波 频率扩展了12倍
            remCodePhase = Tracking.trackResults(channelNr).remCodePhase;
            codePhaseStep = codeFreq / settings.samplingFreq;

            blksize = ceil((12*settings.codeLength-remCodePhase) / codePhaseStep);
            %% pilot channel data

            % Read in the appropriate number of samples to process this
            % interation
            [rawSignal, samplesRead] = fread(fid, ...
                dataAdaptCoeff*blksize, settings.dataType);

            rawSignal = rawSignal';

            if (dataAdaptCoeff==2)
                rawSignal1=rawSignal(1:2:end);
                rawSignal2=rawSignal(2:2:end);
                rawSignal = rawSignal1 - 1i .* rawSignal2;  %transpose vector
            end

            % If did not read in enough samples, then could be out of
            % data - better exit
            if (samplesRead ~= dataAdaptCoeff*blksize)
                disp('Not able to read the specified number of samples  for tracking, exiting!')
                fclose(fid);

                return
            end

            %% Set up all the code phase tracking information -------------------------
            % Update tracking information for code loop

            % Define index into early code vector
            tcode       = (remCodePhase+(earlyLateSpc*12)) : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase+(earlyLateSpc*12));
            tcode2      = ceil(tcode) + 12;
            earlyCode   = pilotChannelCode(tcode2);

            % Define index into late code vector
            tcode       = (remCodePhase-(earlyLateSpc*12)) : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase-(earlyLateSpc*12));
            tcode2      = ceil(tcode) + 12;
            lateCode    = pilotChannelCode(tcode2);

            % Define index into prompt code vector
            tcode       = remCodePhase : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase);
            tcode2      = ceil(tcode) + 12;
            promptCode  = pilotChannelCode(tcode2);
            %Define index into prompt code vector for data channel
            promptCodeData =dataChannelBocCode(tcode2);
            
            remCodePhase_New = (tcode(blksize) + codePhaseStep) - 4092*12;

            % Define index into very early code vector
            tcode       = (remCodePhase+(veryEarlyLateSpc*12)) : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase+(veryEarlyLateSpc*12));
            tcode2      = ceil(tcode) + 12;
            veryEarlyCode   = pilotChannelCode(tcode2);
            % Define index into very late code vector
            tcode       = (remCodePhase-(veryEarlyLateSpc*12)) : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase-(veryEarlyLateSpc*12));
            tcode2      = ceil(tcode) + 12;
            veryLateCode   = pilotChannelCode(tcode2);

            %% Generate the carrier frequency to mix the signal to baseband -----------
            time    = (0:blksize) ./ settings.samplingFreq;
            carrFreq = settings.IF + Tracking.trackResults(channelNr).carrierDoppler;
            remCarrPhase = Tracking.trackResults(channelNr).remCarrPhase;

            % Get the argument to sin/cos functions
            trigarg = ((carrFreq * 2.0 * pi) .* time) + remCarrPhase;
            remCarrPhase_New = rem(trigarg(blksize+1), (2 * pi));

            % Finally compute the signal to mix the collected data to bandband
            carrsig = exp(1i .* trigarg(1:blksize)); %mind this equation
            %% Integrated Doppler

            %% Generate the six standard accumulated values ---------------------------
            % First mix to baseband
            qBasebandSignal = real(carrsig .* rawSignal);
            iBasebandSignal = imag(carrsig .* rawSignal);

            % Now get early, late, and prompt values for each
            I_E = sum(earlyCode  .* iBasebandSignal);
            Q_E = sum(earlyCode  .* qBasebandSignal);
            I_P = sum(promptCode .* iBasebandSignal);
            Q_P = sum(promptCode .* qBasebandSignal);
            I_L = sum(lateCode   .* iBasebandSignal);
            Q_L = sum(lateCode   .* qBasebandSignal);
            I_VE= sum(veryEarlyCode   .* iBasebandSignal);
            Q_VE= sum(veryEarlyCode   .* qBasebandSignal);
            I_VL= sum(veryLateCode   .* iBasebandSignal);
            Q_VL= sum(veryLateCode   .* qBasebandSignal);
            
            data_I_P = sum(promptCodeData .* iBasebandSignal);
            data_Q_P = sum(promptCodeData .* qBasebandSignal);
            %% Locking detection
            %% Find FLL error and update carrier NCO
            % FLL is not stable in this SDR....
            Pdot    = I_P_1*I_P+Q_P_1*Q_P;
            Pcross  = I_P_1*Q_P-Q_P_1*I_P;
%             weError = Pcross*sign(Pdot)/(sqrt(I_P_1^2+Q_P_1^2)*sqrt(I_P^2+Q_P^2))/PDIcarr;
            weError = atan2(Pcross,Pdot)/PDIcarr;
            
            if strcmp(settings.TrackingSettingFLLUsed,'true')
                weNco=FLLLoopFilter(PDIcarr,...
                    settings.FLLa2,settings.FLLBL2,...
                    weError,weOldError_1,weOldError_2,weOldNco_1,weOldNco_2);
            end
       
            %% Find PLL error and update carrier NCO ----------------------------------

            % Implement carrier loop discriminator (phase detector)
            carrError = atan(Q_P / I_P) ; %鉴相误差
            carrNco=loopFilter2(PDIcarr,...
                settings.PLLa3,settings.PLLb3,settings.PLLBL3,...
                carrError,carrOldError_1,carrOldError_2,carrOldNco_1,carrOldNco_2);
            
            %% PLL+FLL carrier frequency controller
            if loopCnt <0
                 a=0;
                 b=1;
            else
                a=1;
                b=0;
            end         
            carrFrqNco=a*0.2*carrNco+b*0.07*weNco;
            carrierDoppler_new = channel(channelNr).acquiredFreq + carrFrqNco;
            %% Find DLL error and update code NCO
%            codeError =(1-earlyLateSpc) *(sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) / ...
%                (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));
           P_early = sqrt(I_VE*I_VE + Q_VE*Q_VE + I_E * I_E + Q_E * Q_E);
           P_late = sqrt(I_VL*I_VL + Q_VL*Q_VL + I_L * I_L + Q_L * Q_L);
           if (P_early + P_late == 0.0)
               codeError =0.0;
           else
               codeError =(P_early - P_late) / (P_early + P_late);
           end
           codeNco=loopFilter1(PDIcarr,...
               settings.DLLa2,settings.DLLBL2,...
               codeError,codeOldError_1,codeOldNco_1);
           codeFrqNco = 2.0*codeNco*12;
           codeDoppler_new = carrierDoppler_new  * 12 * settings.codeFreqBasis/ settings.carrFreqBasis + ...
               codeFrqNco ;

            %% Update tracking information and record various measures to show in postprocessing ----------------------
           
            Tracking.trackResults(channelNr).I_P_1=I_P;
            Tracking.trackResults(channelNr).Q_P_1=Q_P;
            Tracking.trackResults(channelNr).codePhase=1;
            
            Tracking.trackResults(channelNr).carrFreq     = carrFreq;
            Tracking.trackResults(channelNr).remCarrPhase = remCarrPhase_New;
            Tracking.trackResults(channelNr).remCodePhase = remCodePhase_New;
            
            Tracking.trackResults(channelNr).carrierDoppler    = carrierDoppler_new;
            Tracking.trackResults(channelNr).codeDoppler       = codeDoppler_new;
            
            %FLL tracking loop parameters
            Tracking.trackResults(channelNr).trackingweOldNco_2 = Tracking.trackResults(channelNr).trackingweOldNco_1;
            Tracking.trackResults(channelNr).trackingweOldNco_1 = weNco;
            Tracking.trackResults(channelNr).trackingweOldError_2 = Tracking.trackResults(channelNr).trackingweOldError_1;
            Tracking.trackResults(channelNr).trackingweOldError_1 = weError;
            
            %DLL code tracking loop parameters
            Tracking.trackResults(channelNr).trackingcodeOldNco_1   = codeNco;
            Tracking.trackResults(channelNr).trackingcodeOldError_1 = codeError;
            
            %PLL carrier/Costas loop parameters
            Tracking.trackResults(channelNr).trackingcarrOldNco_2   = Tracking.trackResults(channelNr).trackingcarrOldNco_1;
            Tracking.trackResults(channelNr).trackingcarrOldNco_1   = carrNco;
            Tracking.trackResults(channelNr).trackingcarrOldError_2 = Tracking.trackResults(channelNr).trackingcarrOldError_1;
            Tracking.trackResults(channelNr).trackingcarrOldError_1 = carrError;
            %% Record sample number (based on 8bit samples)
            trackResults(channelNr).absoluteSample(loopCnt) = (ftell(fid))/dataAdaptCoeff/settings.Datalength- remCodePhase/codePhaseStep;

            trackResults(channelNr).dllDiscr(loopCnt)       = codeError;
            trackResults(channelNr).dllDiscrFilt(loopCnt)   = codeNco;
            trackResults(channelNr).pllDiscr(loopCnt)       = carrError;
            trackResults(channelNr).pllDiscrFilt(loopCnt)   = carrNco;

            trackResults(channelNr).I_E(loopCnt) = I_E;
            trackResults(channelNr).I_P(loopCnt) = I_P;
            trackResults(channelNr).I_L(loopCnt) = I_L;
            trackResults(channelNr).Q_E(loopCnt) = Q_E;
            trackResults(channelNr).Q_P(loopCnt) = Q_P;
            trackResults(channelNr).Q_L(loopCnt) = Q_L;
            
            trackResults(channelNr).data_I_P(loopCnt) = data_I_P;
            trackResults(channelNr).data_Q_P(loopCnt) = data_Q_P;

            % Evaluate the tracking results status here to ensure the
            % plotTracking to plot the results tracked so far 
            % (In case the tracking update window is closed)
        end % for loopCnt
        trackResults(channelNr).status  = channel(channelNr).status;

        % If we got so far, this means that the tracking was successful
        % Now we only copy status, but it can be update by a lock detector
        % if implemented
      
    end % if a PRN is assigned
end % for channelNr

% Close the waitbar
close(hwb)
