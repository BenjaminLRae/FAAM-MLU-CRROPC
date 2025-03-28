%% Process Data

disp('Processing data for plot...')

% Change length of buffers if number of samples is less than application buffer.

% if(totalValues < appBufferSize)
% 
%     pAppBufferChA.Value(totalValues + 1:end) = [];
%     pAppBufferChB.Value(totalValues + 1:end) = [];
% 
% end

% Retrieve data and convert to milliVolts
% bufferChA = adc2mv(pAppBufferChA.Value, chARangeMv, maxADCValue);
% bufferChB = adc2mv(pAppBufferChB.Value, chBRangeMv, maxADCValue);

% Remove spurious out-of-range data (I don't know what causes this, but
% there appears to be occasional data points with a value of minimum range
for i = 1:length(bufferChA)
    if bufferChA(i) < -0
        bufferChA(i) = 658.132;
    end
end

for i = 1:length(bufferChB)
    if bufferChB(i) < -1999
        bufferChB(i) = -265.084;
    end
end


%% find peaks

% [pks,locs,w,p] = findpeaks(bufferChA,'MinPeakProminence',90);

%% Plot data

figure;
subplot(2,3,1);

% finalFigure = figure('Name','PicoScope 2000 Series Example - Fast Streaming Mode Capture', ...
    % 'NumberTitle','off');

% finalFigureAxes = axes('Parent', finalFigure);
% hold(finalFigureAxes, 'on');
hold on;

% Find the maximum voltage range and add 500mV.
maxYRange = max(chARangeMv, chBRangeMv) + 500;
ylim([(-1 * maxYRange) maxYRange]);

% Calculate time axis in milliseconds, then plot data.
timeLabel = 'Time (ms)';
time = (0:1:(length(bufferChA)) - 1) * double(samplingIntervalMs) * double(numSamplesPerAggregate);
plot(time, bufferChA, time, bufferChB);

title(strcat(string(DateStrings)," ChA: ",detectors{1}," ChB: ",detectors{2}," SInt: ",num2str(samplingIntervalMs*1000000),' ns'));
xlabel(timeLabel);
ylabel('Voltage (mv)');

grid on;
legend('Channel A', 'Channel B','Location','southeast');
hold off;

%% fft
Fs = 1.56e6;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(bufferChA);             % Length of signal
t = (0:L-1)*T;        % Time vector

f = Fs/L*(0:L-1);
fftBufferChA = fft(bufferChA);

subplot(2,3,4);

plot(f,abs(fftBufferChA),"LineWidth",2) 
xlim([0 500]);
ylim([0 8e8]);
title("Single-Sided Amplitude Spectrum of signal")
xlabel("f (Hz)")
ylabel("|fft(bufferChA)|")
grid on;

%% Thresholding to find particle interactions

% Use a threshold of 833.918 mV. Anything above this probably indicates a
% particle interaction. We want to split each of these into a separate
% array so that number of arrays = number of particle interactions. 

threshold_mV_chA = 878;

flag_OverThreshold = false;
flag_UnderThreshold = true;

particleInteractionStartIndex = 0;
particleInteractionEndIndex = 0;

particleInteractionsChA = {};

particleInteractionsCounterChA = 1;

particleInteractionsIndicesChA = [0; 0];

for i = 1:length(bufferChA)
    if flag_UnderThreshold
        % currently we are under the threshold, so in a regime of no
        % particle interaction

        % if this sample is above the threshold, indicates a particle
        % interaction. Therefore set the flag
        if bufferChA(i) > threshold_mV_chA
            flag_OverThreshold = true;
            flag_UnderThreshold = false;

            particleInteractionStartIndex = i;
        end
        
    elseif flag_OverThreshold
        % currently we are over the threshold, so in a regime of current
        % particle interaction. 

        % if this sample is above the threshold, we don't need to do
        % anything - still in the interaction.

        % if this sample is below the threshold, we have completed a
        % particle interaction and need to log the data between start and
        % finish.
        
        if bufferChA(i) > threshold_mV_chA

        elseif bufferChA(i) < threshold_mV_chA
            flag_OverThreshold = false;
            flag_UnderThreshold = true;

            particleInteractionEndIndex = i;

            particleInteractionsChA{particleInteractionsCounterChA} = bufferChA(particleInteractionStartIndex:particleInteractionEndIndex);

            particleInteractionsIndicesChA(1,particleInteractionsCounterChA) = particleInteractionStartIndex;
            particleInteractionsIndicesChA(2,particleInteractionsCounterChA) = particleInteractionEndIndex;
            
            % increase the interactions counter
            particleInteractionsCounterChA = particleInteractionsCounterChA + 1;
        end

    end    

end

%% Thresholding for Channel B

threshold_mV_chB = -283;

flag_OverThreshold = false;
flag_UnderThreshold = true;

particleInteractionStartIndex = 0;
particleInteractionEndIndex = 0;

particleInteractionsChB = {};

particleInteractionsCounterChB = 1;

particleInteractionsIndicesChB = [0; 0];

for i = 1:length(bufferChB)
    if flag_UnderThreshold
        % currently we are under the threshold, so in a regime of no
        % particle interaction

        % if this sample is above the threshold, indicates a particle
        % interaction. Therefore set the flag
        if bufferChB(i) < threshold_mV_chB
            flag_OverThreshold = true;
            flag_UnderThreshold = false;

            particleInteractionStartIndex = i;
        end
        
    elseif flag_OverThreshold
        % currently we are over the threshold, so in a regime of current
        % particle interaction. 

        % if this sample is above the threshold, we don't need to do
        % anything - still in the interaction.

        % if this sample is below the threshold, we have completed a
        % particle interaction and need to log the data between start and
        % finish.
        
        if bufferChB(i) < threshold_mV_chB

        elseif bufferChB(i) > threshold_mV_chB
            flag_OverThreshold = false;
            flag_UnderThreshold = true;

            particleInteractionEndIndex = i;

            particleInteractionsChB{particleInteractionsCounterChB} = bufferChB(particleInteractionStartIndex:particleInteractionEndIndex);

            particleInteractionsIndicesChB(1,particleInteractionsCounterChB) = particleInteractionStartIndex;
            particleInteractionsIndicesChB(2,particleInteractionsCounterChB) = particleInteractionEndIndex;
            
            % increase the interactions counter
            particleInteractionsCounterChB = particleInteractionsCounterChB + 1;
        end

    end    

end

%% Plot each individual particle interaction on same plot
% figure;
% hold on;
% for i = 1:length(particleInteractionsChB)
%     plot(time(particleInteractionsIndicesChB(1,i):particleInteractionsIndicesChB(2,i)),bufferChB(particleInteractionsIndicesChB(1,i):particleInteractionsIndicesChB(2,i)));
% end

%% Measure stats for each particle interaction

particleInteractionsPeakHeightsChA = [];
particleInteractionsTimeWidthsChA = [];
particleInteractionsPeakHeightsChB = [];
particleInteractionsTimeWidthsChB = [];

for i = 1:length(particleInteractionsChA)
    
    % calculate interaction peak height (mV)
    particleInteractionsPeakHeightsChA(i) = max(particleInteractionsChA{i});

    % calculate interaction time (ms)
    particleInteractionsTimeWidthsChA(i) = length(particleInteractionsChA{i}) * samplingIntervalMs;

end

for i = 1:length(particleInteractionsChB)
    
    % calculate interaction peak height (mV)
    particleInteractionsPeakHeightsChB(i) = min(particleInteractionsChB{i});

    % calculate interaction time (ms)
    particleInteractionsTimeWidthsChB(i) = length(particleInteractionsChB{i}) * samplingIntervalMs;

end

%% Plot histograms

subplot(2,3,2);

hChA = histogram(particleInteractionsPeakHeightsChA,64);
xlabel('Interaction peak voltage (mV)');
title('Histogram - peak voltages - Channel A, 64 bins');

subplot(2,3,5);

hChB = histogram(abs(particleInteractionsPeakHeightsChB),64);
xlabel('Interaction peak voltage (mV)');
title('Histogram - peak voltages - Channel B, 64 bins');

%% Calculate normalised concentration for ChA

dNdLogDchA = [];
xVoltageValuesChA = [];

for i = 1:hChA.NumBins

    xVoltageValuesChA(i) = hChA.BinEdges(i) + ((hChA.BinEdges(i+1) - hChA.BinEdges(i))/2);
    
    dN = hChA.Values(i);
    dlogD = log(hChA.BinEdges(i+1)) - log(hChA.BinEdges(i));
    dNdLogDchA(i) = dN / dlogD;

end

subplot(2,3,3);

plot(xVoltageValuesChA,dNdLogDchA,"square-");
title("Voltage distribution of particles - Ch A, 64 Bins, 10s");
xlabel('Interaction peak voltage (mV)');
ylabel('dN/dlogmV (cm^-^3)')
xscale log;
% yscale log;
grid on;

%% Calculate normalised concentration for ChB

dNdLogDchB = [];
xVoltageValuesChB = [];

for i = 1:hChB.NumBins

    xVoltageValuesChB(i) = hChB.BinEdges(i) + ((hChB.BinEdges(i+1) - hChB.BinEdges(i))/2);
    
    dN = hChB.Values(i);
    dlogD = log(hChB.BinEdges(i+1)) - log(hChB.BinEdges(i));
    dNdLogDchB(i) = dN / dlogD;

end

subplot(2,3,6);

plot(abs(xVoltageValuesChB),abs(dNdLogDchB),"square-");
title("Voltage distribution of particles - Ch B, 64 Bins, 10s");
xlabel('Interaction peak voltage (mV)');
ylabel('dN/dlogmV (cm^-^3)')
xscale log;
% yscale log;
grid on;


%% Convert mV values to light power values using detector sensitivity

sensitivityChA = 1.25e+8;
particleInteractionsPeakHeightsChA_W = [];

for i = 1:length(particleInteractionsPeakHeightsChA)

    particleInteractionsPeakHeightsChA_W(i) = (particleInteractionsPeakHeightsChA(i) * 1e-3) / sensitivityChA;

end