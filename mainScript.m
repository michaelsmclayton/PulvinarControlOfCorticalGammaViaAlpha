clear all; close all; clc;
addpath(genpath('/Users/michaelclayton/Documents/Neuroscience/ComputationalModels/chronux_2_12'))

%% Simulation parameters

% Simulation length
numberOfTrials = 5;
fullSimulationLength = 2000;
windowToRemove = 50; % The first few milliseconds are quite unsettled
simulationLength = fullSimulationLength + windowToRemove;
timePoints = linspace(1, simulationLength, simulationLength);

% Oscillator parameters
phaseShift = 0;
alphaAmplitude = 3;

% Plot parameters
fontSize = 14;

% Initialise data stores
allAreaFirings = cell(numberOfTrials, 1);

%% Simulation loop
for trial = 1:numberOfTrials
    
    disp(['Simulation trial = ' num2str(trial)])

    % Create NeuralAreas
    area1 = NeuralArea(simulationLength);
    area2 = NeuralArea(simulationLength);

    % Create Oscillators

    oscillator1 = Oscillator(0, alphaAmplitude, simulationLength);
    oscillator2 = Oscillator(phaseShift, alphaAmplitude, simulationLength);
    
    % Single trial
    for t = timePoints

        % Update oscillators
        oscillator1.update(t)
        oscillator2.update(t)

        % Update neural areas
        area1.update(t, oscillator1.currentVoltage, area2)
        area2.update(t, oscillator2.currentVoltage, false)
    end

    % Remove opening window from sample
    area1.firings(area1.firings(:,1)<=windowToRemove,:) = [];
    area2.firings(area2.firings(:,1)<=windowToRemove,:) = [];
    area1.firings(:,1) = area1.firings(:,1) - windowToRemove;
    area2.firings(:,1) = area2.firings(:,1) - windowToRemove;
    
    allAreaFirings{trial} = {area1.firings, area2.firings};
    
end


%% Calculate and Spike Time Histogram (STM) and Power
STM = {};
areaToAnalyse = 1;

% Frequency spectrum parameters
params.Fs = 1000;
params.fpass = 4:100;
params.tapers = [5 10];
params.trialave = 1;
params.err = [2, .05];

% Spectrogram parameters
window = [.3 0.02];
specParams = params;
specParams.fpass = 5:100;
specParams.pad = -1;
specParams.tapers= [4 window(1) 1];

% Wavelet parameters
scales = logspace(log10(150),.8);
F = scal2frq(scales,'cmor1-1',1/1000);

% Initiaise storage
powerStore = []; %zeros(numberOfTrials, length(params.fpass)+1);
spectrogramStore = [];

% Loop
for trial = 1:numberOfTrials

    % Calculate STM
    data = allAreaFirings{trial}{areaToAnalyse};
    for t = 1:simulationLength
        currentIndices = find(data(:,1)==t);
        if not(isempty(currentIndices))
            excitatorySpikes = find(data(currentIndices,2)<=400);
            inhibitorySpikes = find(data(currentIndices,2)>400);
            STM{1}(t) = length(excitatorySpikes);
            STM{2}(t) = -length(inhibitorySpikes);
        end
    end
    
    % Calculate power
    [powerStore(trial,:), fPower, Serr]= mtspectrumc(STM{1},params);
    
    % Calculate current spectrogram
    currentCoefficients = cwt(STM{1}, scales,'cmor1-1','plot');
    spectrogramStore(trial,:,:) = wscalogram('image', squeeze(currentCoefficients), 'scales', F);
    
    % Normalising and smoothing
    smoothing = 10;
    for f = 1:size(spectrogramStore(trial,:,:),2)
        spectrogramStore(trial,f,:) = spectrogramStore(trial,f,:) / max(spectrogramStore(trial,f,:));
        spectrogramStore(trial,f,:) = smooth(spectrogramStore(trial,f,:), smoothing);
    end
    
    % [spectrogramStore(trial,:,:), times, fSpec, err] = mtspecgramc(STM{1}', window, specParams);

end


%% Plot spectrograms
for i = 1:numberOfTrials
    figure(200+i)
    data = squeeze(spectrogramStore(i,:,:));
    data = flipud(data);
    imagesc(data)
    caxis([0 1]);
    colorbar;
    colormap jet
    ax = gca;
    ax.YTick = 1:5:length(F);
    ax.YTickLabel = fliplr(F(ax.YTick));
    xlim([1 1000])
end

% % % Plot power
% % figure(5)
% % plot(fPower, mean(powerStore))
% 
% % Plot wavelet
% for i = 1:numberOfTrials
%     figure(200+i)
%     wscalogram('image', squeeze(spectrogramStore(i,:,:)), 'scales', F);
%     ylim([2 50])
%     title(' ')
%     ylabel('Frequency (Hz)')
%     xlabel('Time (ms)')
%     colormap default
% end


%% Plot rastergram of excitatory (red) and inhibitory (blue) neurons in one population during an interval of 1000 ms.
hFig = figure(1); hold on;
set(hFig, 'Position', [10 10 600 500]) 
excitatorySpikes1 = find(area1.firings(:,2)<=area1.n_regularSpiking);
inhibitorySpikes1 = find(area1.firings(:,2)>area1.n_regularSpiking);
excitatorySpikes2 = find(area2.firings(:,2)<=area2.n_regularSpiking);
inhibitorySpikes2 = find(area2.firings(:,2)>area2.n_regularSpiking);
plot(area1.firings(excitatorySpikes1,1), area1.firings(excitatorySpikes1,2), '.', 'Color', 'r');
plot(area1.firings(inhibitorySpikes1,1), area1.firings(inhibitorySpikes1,2), '.', 'Color', [0 .2 .7]);
ylabel('Neuron #')
xlabel('Time (ms)') 
set(gca,'FontSize', fontSize)
% subplot(2,1,2); hold on;
% plot(area2.firings(excitatorySpikes2,1), area2.firings(excitatorySpikes2,2), '.', 'Color', 'r');
% plot(area2.firings(inhibitorySpikes2,1), area2.firings(inhibitorySpikes2,2), '.', 'Color', [0 .2 .7]);
% set(gca,'FontSize', fontSize)


%%  Plot Spike Time Histogram
red = [1 .2 0]; blue = [0 .2 1];
hFig = figure(2); hold on;
set(hFig, 'Position', [10 10 600 500])
smoothing = 10;
% subplot(2,1,a); hold on;
area(smooth(STM{1}, smoothing), 'FaceColor', red, 'LineStyle', 'none');
area(smooth(STM{2}, smoothing), 'FaceColor', blue, 'LineStyle', 'none');
plot([0 1000], [0 0], 'Color', [0 0 0], 'LineWidth', 1)
ylabel('# Neurons firing')
xlabel('Time (ms)') 
set(gca,'FontSize', fontSize)
ylim([-7 10])


% Plot oscillators
% figure(1); hold on
% plot(oscillator1.timeseries)
% plot(oscillator2.timeseries)





