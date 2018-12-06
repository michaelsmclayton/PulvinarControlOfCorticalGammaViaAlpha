clear all; close all; clc;

%% Stimulation parameters
simulationLength = 1000; % simulation of 1000 ms
timePoints = linspace(1, simulationLength, simulationLength);

% Create NeuralAreas
area1 = NeuralArea();
area2 = NeuralArea();

% Create Oscillators
phaseShift = 30;
oscillator1 = Oscillator(0);
oscillator2 = Oscillator(phaseShift);

% Simulation loop
for t = timePoints
    
    % Update oscillators
    oscillator1.update(t)
    oscillator2.update(t)
    
    % Update neural areas
    area1.update(t, oscillator1.alphaAmplitude, area2)
    area2.update(t, oscillator2.alphaAmplitude, area1)
    
end


%% Plot oscillators
% figure(1); hold on
% plot(oscillator1.timeseries)
% plot(oscillator2.timeseries)

%% Plot graph showing all spikes during simulation window
figure(100); hold on;
excitatorySpikes = find(area1.firings(:,2)<=area1.n_regularSpiking);
inhibitorySpikes = find(area1.firings(:,2)>area1.n_regularSpiking);
plot(area1.firings(excitatorySpikes,1), area1.firings(excitatorySpikes,2), '.', 'Color', 'r');
plot(area1.firings(inhibitorySpikes,1), area1.firings(inhibitorySpikes,2), '.', 'Color', [0 .2 .7]);


% %%
% 
% %% asdas
% % Time specifications:
% Fs = 1000; % samples per second
% dt = 1/Fs; % seconds per sample
% t = linspace(, 1000, 1000);
% 
% phaseDifference = .2;
% 
% % Sine wave:
% frequency = 10; % Hertz
% % x = cos(2*pi*frequency*t);
% 





