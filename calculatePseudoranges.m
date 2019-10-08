function [pseudoranges] = calculatePseudoranges(trackResults, ...
                                                msOfTheSignal, ...
                                                channelList, settings)
%calculatePseudoranges finds relative pseudoranges for all satellites
%listed in CHANNELLIST at the specified millisecond of the processed
%signal. The pseudoranges contain unknown receiver clock offset. It can be
%found by the least squares position search procedure. 
%
%[pseudoranges] = calculatePseudoranges(trackResults, msOfTheSignal, ...
%                                       channelList, settings)
%
%   Inputs:
%       trackResults    - output from the tracking function
%       msOfTheSignal   - pseudorange measurement point (millisecond) in
%                       the trackResults structure
%       channelList     - list of channels to be processed
%       settings        - receiver settings
%
%   Outputs:
%       pseudoranges    - relative pseudoranges to the satellites. 

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis
% Written by Darius Plausinaitis
% Based on Peter Rinder and Nicolaj Bertelsen
%--------------------------------------------------------------------------
%
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

% CVS record:
% $Id: calculatePseudoranges.m,v 1.1.2.18 2006/08/09 17:20:11 dpl Exp $

%--- Set initial travel time to infinity ----------------------------------
% Later in the code a shortest pseudorange will be selected. Therefore
% pseudoranges from non-tracking channels must be the longest - e.g.
% infinite. 
travelTime = inf(1, settings.numberOfChannels);

% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
                        (settings.codeFreqBasis / settings.codeLength));

%--- For all channels in the list ... 
for channelNr = 1:numel(channelList) % numel required becase sometime row, sometimes column
    %--- Compute the travel times -----------------------------------------    
    travelTime(channelList(channelNr)) = ...
        4 * trackResults(channelList(channelNr)).absoluteSample((msOfTheSignal(channelList(channelNr)))) / samplesPerCode;
end

%--- Truncate the travelTime and compute pseudoranges ---------------------
% http://www.insidegnss.com/node/2898
% This procedure yields relative arrival times between satellites, but not the absolute pseudoranges. To determine the set of pseudoranges for the first time, the channel with the earliest arriving subframe is assumed as reference and a minimum travel time is assigned to the reference channel based on the known orbits of the satellites and typical user altitudes (e.g., its value, for GPS, is in a range between 65 and 85 milliseconds).
% All other pseudoranges are then derived with respect to the reference channel by adding the relative-arrival times.
% TODO: is the settings offset updated once we know the closest satellite true range? or are we always carrying error in diff between settings.startOffset and truth?
minimum         = floor(min(travelTime));
travelTime      = travelTime - minimum + settings.startOffset;

%--- Convert travel time to a distance ------------------------------------
% The speed of light must be converted from meters per second to meters
% per millisecond. 
pseudoranges    = travelTime * (settings.c / 1000);
