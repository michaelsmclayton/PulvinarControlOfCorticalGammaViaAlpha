clear all; close all; clc;

%% Simulation
simulationLength = 1000; % simulation of 1000 ms
area1 = NeuralArea(simulationLength);
timeStep = 1;
timePoint = 1;
for t = 1:timeStep:simulationLength
    area1.update(t, timePoint)
    timePoint = timePoint + 1;
end

%% Plot graph showing all spikes during simulation window
figure(1); hold on;
excitatorySpikes = find(area1.firings(:,2)<=area1.n_regularSpiking);
inhibitorySpikes = find(area1.firings(:,2)>area1.n_regularSpiking);
plot(area1.firings(excitatorySpikes,1), area1.firings(excitatorySpikes,2), '.', 'Color', 'r');
plot(area1.firings(inhibitorySpikes,1), area1.firings(inhibitorySpikes,2), '.', 'Color', [0 .2 .7]);


