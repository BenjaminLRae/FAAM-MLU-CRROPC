%% PicoScope 2000 Series Instrument Driver Oscilloscope Block Data Capture Example 
% This is a modified version of the machine generated representation of an
% instrument control session using a device object. The instrument control
% session comprises all the steps you are likely to take when communicating
% with your instrument.
% 
% These steps are:
%       
% # Create a device object   
% # Connect to the instrument 
% # Configure properties 
% # Invoke functions 
% # Disconnect from the instrument 
%  
% To run the instrument control session, type the name of the file,
% PS2000_ID_Block_Example, at the MATLAB command prompt.
% 
% The file, PS2000_ID_BLOCK_EXAMPLE.M must be on your MATLAB PATH. For
% additional information on setting your MATLAB PATH, type 'help addpath'
% at the MATLAB command prompt.
%
% *Example:*
%     PS2000_ID_Block_Example;
%
% *Description:*
%     Demonstrates how to call functions in order to collect a block of 
%     data from a PicoScope 2000 Series Oscilloscope.
%
% *See also:* <matlab:doc('icdevice') |icdevice|> | <matlab:doc('instrument/invoke') |invoke|>
%
% *Copyright:* © 2013 - 2017 Pico Technology Ltd. All rights reserved.

%% Suggested Input Test Signals
% This example was published using the following test signals:
%
% * Channel A: 4Vpp, 5kHz sine wave
% * Channel B: 2Vpp, 5kHz square wave

%% Clear Command Window and Close any Figures

clc;
close all;

%% Load Configuration Information

PS2000Config;

%% Device Connection

% Create a device object. 
ps2000DeviceObj = icdevice('picotech_ps2000_generic.mdd');

% Connect device object to hardware.
connect(ps2000DeviceObj);

%% Obtain Device Groups
% Obtain references to device groups to access their respective properties
% and functions.

% Block specific properties and functions are located in the Instrument
% Driver's Block group.

blockGroupObj = get(ps2000DeviceObj, 'Block');
blockGroupObj = blockGroupObj(1);

%% Device Configuration
% Default channel setup, sampling interval and number of samples for the
% device are used.
%
% Refer to the 'Configure Device' section in the
% <PS2000_ID_Block_SimpleTrig_Example.html PicoScope 2000 Series Instrument Driver Oscilloscope Block Data Capture with Simple Trigger Example> 
% for an example of how to set the channels, sampling interval and number
% of samples.

%%
% Set channels:

% Channel     : 0 (ps2000Enuminfo.enPS2000Channel.PS2000_CHANNEL_A)
% Enabled     : 1 (PicoConstants.TRUE)
% DC          : 1 (DC Coupling)
% Range       : 6 (ps2000Enuminfo.enPS2000Range.PS2000_1V)
[status.setChA] = invoke(ps2000DeviceObj, 'ps2000SetChannel', 0, 1, 1, 8);    

% Channel     : 1 (ps2000Enuminfo.enPS2000Channel.PS2000_CHANNEL_B)
% Enabled     : 1 (PicoConstants.TRUE)
% DC          : 1 (DC Coupling)
% Range       : 7 (ps2000Enuminfo.enPS2000Range.PS2000_2V)
[status.setChB] = invoke(ps2000DeviceObj, 'ps2000SetChannel', 1, 1, 1, 8);  

%%
% Set sampling interval and number of samples to collect:

% Set sampling interval - the actual sampling interval used by the driver
% will be returned.

% timeIntervalUs : 100 (microseconds)
[samplingIntervalUs, maxBlockSamples] = invoke(blockGroupObj, 'setBlockIntervalUs', 100);

% Confirm the timebase index selected.
timebaseIndex = get(blockGroupObj, 'timebase'); 
disp(strcat('Timebase Index selected: ',num2str(timebaseIndex)));

totalSamplingTime = 0.1;

numStreamingSamples = (1000 * totalSamplingTime) / samplingIntervalUs;

% Set the number of samples to collect.
set(ps2000DeviceObj, 'numberOfSamples', 3968);

%% Data Collection
% Capture a block of data on Channels A and B together with times. Data for
% channels is returned in millivolts.
disp('Collecting block of data...');

% Execute device object function(s).
[bufferTimes, bufferChA, bufferChB, numDataValues, timeIndisposedMs] = invoke(blockGroupObj, 'getBlockData');

disp('Data collection complete.');

%% Stop the Device
% Additional blocks can be captured prior to stopping the device.

stopStatus = invoke(ps2000DeviceObj, 'ps2000Stop');

%% Process Data
% Process data as required. In this example the data is displayed in a
% figure.

disp('Plotting data...')

% Find the time units used by the driver.
timesUnits = timeunits(get(blockGroupObj, 'timeUnits'));

% Append to string.
timeLabel = strcat('Time (', timesUnits, ')');

% Plot the data.

figure1 = figure('Name','PicoScope 2000 Series Example - Block Mode Capture', ...
    'NumberTitle', 'off');

plot(bufferTimes, bufferChA, bufferTimes, bufferChB);
title('Block Data Acquisition');
xlabel(timeLabel);
ylabel('Voltage (mv)');
legend('Channel A', 'Channel B');
grid on;

%% Disconnect
% Disconnect device object from hardware.
disconnect(ps2000DeviceObj);
delete(ps2000DeviceObj);
