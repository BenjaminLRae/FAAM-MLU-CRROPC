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
% Fs = 1.56e6;            % Sampling frequency                    
% T = 1/Fs;             % Sampling period       
% L = length(bufferChA);             % Length of signal
% t = (0:L-1)*T;        % Time vector
% 
% f = Fs/L*(0:L-1);
% fftBufferChA = fft(bufferChA);
% 
% subplot(2,3,4);
% 
% plot(f,abs(fftBufferChA),"LineWidth",2) 
% xlim([0 500]);
% ylim([0 8e8]);
% title("Single-Sided Amplitude Spectrum of signal")
% xlabel("f (Hz)")
% ylabel("|fft(bufferChA)|")
% grid on;

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

threshold_mV_chB = -266;

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

%% Convert mV values to light power values using detector sensitivity

sensitivityChA = 1.25e+8;
particleInteractionsPeakHeightsChA_W = [];

for i = 1:length(particleInteractionsPeakHeightsChA)

    particleInteractionsPeakHeightsChA_W(i) = (particleInteractionsPeakHeightsChA(i) * 1e-3) / sensitivityChA;

end

%% Convert light power values to particle size using estimated scattering data

% Import data from text file
% Script for importing data from the following text file:
%
%    filename: G:\My Drive\PCASP\Optics simulation\opc_double_parabola_scattering_values_larger_range3.csv
%
% Auto-generated by MATLAB on 26-Mar-2025 11:20:16

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [7, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Diameter", "ScatteringCrossSection_um"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
scatteringDataTable = readtable("G:\My Drive\PCASP\Optics simulation\opc_double_parabola_scattering_values_larger_range3.csv", opts);

% Clear temporary variables
clear opts

% add column for scattered power
scatteringDataTable.ScatteringCrossSection_m = scatteringDataTable.ScatteringCrossSection_um * 1e-12;
scatteringDataTable.HalfCrossSection_m = scatteringDataTable.ScatteringCrossSection_m / 2; % due to the mirror optics configuration
scatteringDataTable.ScatteredPower = zeros(height(scatteringDataTable),1);

laserPower = 1; % W
laserDiameter = 1e-3; % m
laserBeamArea = pi() * (laserDiameter/2)^2; % m^2
laserIrradiance = laserPower / laserBeamArea; % W/m^2

for i = 1:height(scatteringDataTable)
    
    scatteringDataTable.ScatteredPower(i) = laserIrradiance * scatteringDataTable.HalfCrossSection_m(i);

end

%% interpolation/finding corresponding particle size from scattering data

%vq = interp1(scatteringDataTable.ScatteredPower,scatteringDataTable.Diameter,3e-5,'nearest');

% figure;
% plot(scatteringDataTable.ScatteredPower,scatteringDataTable.Diameter);
% grid on;
% hold on;
% xline(3e-5);
% uniqueScatterPower = unique(scatteringDataTable.ScatteredPower);
% 
% %plot(unique(scatteringDataTable.ScatteredPower),scatteringDataTable.Diameter);

particleInteractionsEstimatedDiameterChA = [];

for i = 1:length(particleInteractionsPeakHeightsChA_W)
    
    flag_foundMatch = false;
    flag_matchUnderRange = false;
    flag_matchOverRange = false;
    returnDiameter = 0;

    for j = 1:height(scatteringDataTable.ScatteredPower)
        
        if particleInteractionsPeakHeightsChA_W(i) < scatteringDataTable.ScatteredPower(j)
            
            % we've found the first instance where the scattered power is
            % higher than the interaction power, so the previous point is
            % the closest fit for particle diameter without going over.

            flag_foundMatch = true;
            
            if j == 1
                
                % scattered power from particle is lower than any simulated
                % values, so it is smaller than 0.01 um.

                flag_matchUnderRange = true;
                returnDiameter = scatteringDataTable.Diameter(j);
                break;

            else

                % scattered power is within simulated values

                returnDiameter = scatteringDataTable.Diameter(j-1);
                break;
            end                        
        end        
    end

    if flag_foundMatch == false
            
            % we didn't find a match, so our scattered power is higher than
            % any simulated values, i.e. it's larger than 15 um.

            flag_matchOverRange = true;
            returnDiameter = scatteringDataTable.Diameter(height(scatteringDataTable.ScatteredPower));
    end

    particleInteractionsEstimatedDiameterChA(i) = returnDiameter;

end

%% Plot histograms

subplot(2,3,2);

hChA = histogram(particleInteractionsPeakHeightsChA_W,64);
xlabel('Interaction peak light power level (W)');
title('Histogram - peak power - Channel A, 64 bins');
grid on;

subplot(2,3,5);

hChA_D = histogram(particleInteractionsEstimatedDiameterChA);
xlabel('Estimated particle diameter (um)');
title(strcat("Histogram - est. diameter - Channel A, ",num2str(hChA_D.NumBins)," bins"));
grid on;

%% Calculate normalised concentration for ChA light power

dNdLogDchA = [];
xPowerValuesChA = [];

for i = 1:hChA.NumBins

    xPowerValuesChA(i) = hChA.BinEdges(i) + ((hChA.BinEdges(i+1) - hChA.BinEdges(i))/2);
    
    dN = hChA.Values(i);
    dlogD = log(hChA.BinEdges(i+1)) - log(hChA.BinEdges(i));
    dNdLogDchA(i) = dN / dlogD;

end

subplot(2,3,3);

plot(xPowerValuesChA,dNdLogDchA,"square-");
title("Light power distribution of particles - Ch A, 64 Bins, 10s");
xlabel('Interaction peak light power (W)');
ylabel('dN/dlogW (cm^-^3)')
xscale log;
% yscale log;
grid on;

%% Calculate normalised concentration for ChA particle diameter

dNdLogDchA_D = [];
xDiameterValuesChA_D = [];

for i = 1:hChA_D.NumBins

    xDiameterValuesChA_D(i) = hChA_D.BinEdges(i) + ((hChA_D.BinEdges(i+1) - hChA_D.BinEdges(i))/2);
    
    dN = hChA_D.Values(i);
    dlogD = log(hChA_D.BinEdges(i+1)) - log(hChA_D.BinEdges(i));
    dNdLogDchA_D(i) = dN / dlogD;

end

subplot(2,3,6);

plot(abs(xDiameterValuesChA_D),abs(dNdLogDchA_D),"square-");
title(strcat("Norm. concentration ChA, ",num2str(hChA_D.NumBins)," bins"));
xlabel('Estimated particle diameter (um)');
ylabel('dN/dlogD (cm^-^3)')
xscale log;
% yscale log;
grid on;

%% Plot particle diameter range

subplot(2,3,4);

x_points = [0, 0, max(scatteringDataTable.Diameter), max(scatteringDataTable.Diameter)];  
y_points = [min(particleInteractionsPeakHeightsChA_W), max(particleInteractionsPeakHeightsChA_W), max(particleInteractionsPeakHeightsChA_W), min(particleInteractionsPeakHeightsChA_W)];
color = [0, 0, 1];

hold on;
a = fill(x_points, y_points, color);
a.FaceAlpha = 0.1;

plot(scatteringDataTable.Diameter,scatteringDataTable.ScatteredPower);
grid on;
xlabel('Particle diameter (um)');
ylabel('Scattered power (W)');

%hold on;
%ha = area([0 max(scatteringDataTable.Diameter)], [min(particleInteractionsPeakHeightsChA_W) max(particleInteractionsPeakHeightsChA_W)]);

%yline(min(particleInteractionsPeakHeightsChA_W));
%yline(max(particleInteractionsPeakHeightsChA_W));


