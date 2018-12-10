clear all; close all; clc;

%% Add necessary paths
addpath(genpath('./functions')) % for helper function
addpath(genpath('./classes')) % for class defintions
addpath(genpath('../chronux_2_12')) % for Chronux library
fontSize = 14; % For subsequent figures

%% Run simulation

% Define simulation parameters
simParams.numberOfTrials = 10;
fullSimulationLength = 2000;
simParams.windowToRemove = 50; % The first few milliseconds are quite unsettled
simParams.simulationLength = fullSimulationLength + simParams.windowToRemove;

% Define oscillator parameters
simParams.phaseShift = 0;
simParams.alphaAmplitude = 2;

% Run model
[area1, area2, oscillator1, oscillator2, allSpikes] = runModel(simParams);

% Calculate and plot sample spike time histogram (STM)
[STM] = calculateSTM(allSpikes, simParams); % shape {trials}(struct: excitatory vs. inhibitory)


%% Calculate power and spectrogram data

% Parameters for power analysis
powerParams.Fs = 1000;
powerParams.fpass = 4:100;
powerParams.tapers = [5 10];
powerParams.trialave = 1;
powerParams.err = [2, .05];
powerParams.scales = logspace(1, log10(150)); % frequencies for spectrogram analysis

% Get data
[excitatoryPowerStore, inhibitoryPowerStore, fPower, ...
    excitatorySpectrogramStore, inhibitorySpectrogramStore, F] = calculatePower(STM, powerParams);


%% Plot power and spectrograms

% % Plot power
% hFig = figure(3);
% set(hFig, 'Position', [10 10 600 500])
% plot(fPower, squeeze(mean(excitatoryPowerStore(:,1,:),1)), 'LineWidth', 2)

% Plot sample spectrogram
hFig = figure(4); hold on;
set(hFig, 'Position', [10 10 600 500])
sampleSpectrogramData = squeeze(excitatorySpectrogramStore(1,7,:,:));
imagesc(flipud(sampleSpectrogramData))
caxis([0 1]);
colorbar;
colormap jet
ax = gca;
ax.YTick = 1:5:length(F);
ax.YTickLabel = fliplr(round(F(ax.YTick)));
ylim([1 50])
xlim([1 1000])
ylabel('Frequency (Hz')
xlabel('Time (ms')
set(gca,'FontSize', fontSize)


%% Analyse gamma power by alpha phase bin

hFig = figure(5); hold on;
set(hFig, 'Position', [10 10 600 500])

% Parameters for gamma power analysis
gammaBand = [20 50]; % Gamma band
gammaIndices = and(F>=gammaBand(1), F<=gammaBand(2));

% Parameters for alpha phase analysis
numberOfBins = 20;
angleTimeSeries = mod(1:1000, 100) / 100;

% Loop over excitatory and inhibitory data
currentArea = 1;
for currentData = 1:2
    
    % Select data
    if currentData == 1
        dataToAnalyse = squeeze(excitatorySpectrogramStore(currentArea,:,:,:));
    else
        dataToAnalyse = squeeze(inhibitorySpectrogramStore(currentArea,:,:,:));
    end

    % Loop through bins
    powerByBin = zeros(numberOfBins-1,1);
    errByBin = zeros(numberOfBins-1,1);
    currentRange = [0 1/numberOfBins];
    for b = 1:(numberOfBins-1)
        currentIndices = find(and(angleTimeSeries>currentRange(1), angleTimeSeries<currentRange(2)));
        currentData = mean(mean(dataToAnalyse(:, gammaIndices', currentIndices'),1),2);
        powerByBin(b) = mean(currentData);
        errByBin(b) = std(currentData)/sqrt(length(currentData)); 
        currentRange = currentRange + 1/numberOfBins;
    end

    % Shift to center on 0 degrees
    powerByBin = vertcat(powerByBin(7:end), powerByBin(1:6));
    errByBin = vertcat(errByBin(7:end), errByBin(1:6));

    % Plot results
    errorbar(powerByBin, errByBin, 'LineWidth', 2);  
end

% Add alpha phase angle information
angles = linspace(-180, 180, 10);
ax = gca;
ax.XTickLabel = angles;


%% Plot rastergram of excitatory (red) and inhibitory (blue) neurons in one population during an interval of 1000 ms.
hFig = figure(1); hold on;
set(hFig, 'Position', [10 10 600 500]) 
excitatorySpikes1 = find(area1.firings(:,2)<=area1.n_regularSpiking);
inhibitorySpikes1 = find(area1.firings(:,2)>area1.n_regularSpiking);
plot(area1.firings(excitatorySpikes1,1), area1.firings(excitatorySpikes1,2), '.', 'Color', 'r');
plot(area1.firings(inhibitorySpikes1,1), area1.firings(inhibitorySpikes1,2), '.', 'Color', [0 .2 .7]);
ylabel('Neuron #')
xlabel('Time (ms)') 
set(gca,'FontSize', fontSize)


%% Plot a sample STM
hFig = figure(2); hold on;
set(hFig, 'Position', [10 10 600 500])
red = [1 .2 0]; blue = [0 .2 1]; smoothing = 10;
excitatoryData = smooth(STM{1}.excitatory(1,:), smoothing);
inhibitoryData = smooth(STM{1}.inhibitory(1,:), smoothing);
area(excitatoryData, 'FaceColor', red, 'LineStyle', 'none');
area(inhibitoryData, 'FaceColor', blue, 'LineStyle', 'none');
plot([0 1000], [0 0], 'Color', [0 0 0], 'LineWidth', 1)
ylabel('# Neurons firing')
xlabel('Time (ms)') 
set(gca,'FontSize', fontSize)
ylim([-7 10])

