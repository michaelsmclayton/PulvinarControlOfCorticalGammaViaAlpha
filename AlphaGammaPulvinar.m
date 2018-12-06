clear all; close all; clc;
% Created by Michael Clayton 6th December 2018
% Tested on MATLAB R2014b

% A MATLAB implementation of 'Top-down control of cortical gamma-band communication
% via pulvinar induced phase shifts in the alpha rhythm [Quax, Jensen, & Tiesinga(2017;
% Computational Biology)

%% Model parameters

% Numbers of each neuron type
n_regularSpiking = 400;
n_fastSpiking = 75;
n_lowThreshSpiking = 25;
totalNeurons = n_regularSpiking + n_fastSpiking + n_lowThreshSpiking;

% Neuron parameters
params.regularSpiking.a = .02; % for Regular spiking neurons
params.regularSpiking.b = .2;
params.regularSpiking.c = -65;
params.regularSpiking.d = 8;
params.fastSpiking.a = .1; % for Fast spiking neurons
params.fastSpiking.b = .2;
params.fastSpiking.c = -65;
params.fastSpiking.d = 2;
params.lowThreshSpiking.a = .02; % for Low-threshold spiking neurons
params.lowThreshSpiking.b = .25;
params.lowThreshSpiking.c = -65;
params.lowThreshSpiking.d = 2;

% Connectivity parameters
connectivity.rs.rs = .0375; % from Regular spiking to...
connectivity.rs.fs = -.25;
connectivity.rs.lts = -.3;
connectivity.fs.rs = .125; % from Fast spiking to...
connectivity.fs.fs = -.15;
connectivity.fs.lts = -.1;
connectivity.lts.rs = .125; % from Low-threshold spiking to...
connectivity.lts.fs = -.1;
connectivity.lts.lts = 0;

% Noise parameters
noise.rs.mean = 4; % for Regular spiking neurons
noise.rs.std = 6;
noise.fs.mean = 5; % for Fast spiking neurons
noise.fs.std = 4;
noise.lts.mean = 4; % for Low-threshold spiking neurons
noise.lts.std = 4;


%% Create the network

% Create arrays of parameters for each neuron
A = vertcat(params.regularSpiking.a * ones(n_regularSpiking,1), ...
            params.fastSpiking.a * ones(n_fastSpiking,1), ...
            params.lowThreshSpiking.a * ones(n_lowThreshSpiking,1));
B = vertcat(params.regularSpiking.b * ones(n_regularSpiking,1), ...
            params.fastSpiking.b * ones(n_fastSpiking,1), ...
            params.lowThreshSpiking.b * ones(n_lowThreshSpiking,1));
C = vertcat(params.regularSpiking.c * ones(n_regularSpiking,1), ...
            params.fastSpiking.c * ones(n_fastSpiking,1), ...
            params.lowThreshSpiking.c * ones(n_lowThreshSpiking,1));
D = vertcat(params.regularSpiking.d * ones(n_regularSpiking,1), ...
            params.fastSpiking.d * ones(n_fastSpiking,1), ...
            params.lowThreshSpiking.d * ones(n_lowThreshSpiking,1));
        
% Create synaptic connections
synapses.rs.rs = connectivity.rs.rs * ones(n_regularSpiking, n_regularSpiking);
synapses.rs.fs = connectivity.rs.fs * ones(n_regularSpiking, n_fastSpiking);
synapses.rs.lts = connectivity.rs.lts * ones(n_regularSpiking, n_lowThreshSpiking);
synapses.fs.rs = connectivity.fs.rs * ones(n_fastSpiking, n_regularSpiking);
synapses.fs.fs = connectivity.fs.fs * ones(n_fastSpiking, n_fastSpiking);
synapses.fs.lts = connectivity.fs.lts * ones(n_fastSpiking, n_lowThreshSpiking);
synapses.lts.rs = connectivity.lts.rs * ones(n_lowThreshSpiking, n_regularSpiking);
synapses.lts.fs = connectivity.lts.fs * ones(n_lowThreshSpiking, n_fastSpiking);
synapses.lts.lts = connectivity.lts.lts * ones(n_lowThreshSpiking, n_lowThreshSpiking);
allSynapses = vertcat(horzcat(synapses.rs.rs, synapses.rs.fs, synapses.rs.lts), ...
                      horzcat(synapses.fs.rs, synapses.fs.fs, synapses.fs.lts), ...
                      horzcat(synapses.lts.rs, synapses.lts.fs, synapses.lts.lts));
                  
% Set, for all neurons, initial value of V to -65 (membrane potential of the simulated neuron)
v = -65*ones(totalNeurons,1);

% Set 'u' values as initial 'v' values, multipled by 'b' (sensitivity to 'u')
u = B.*v; % Initial values of u


%% Simulation
simulationLength = 1000; % simulation of 1000 ms
v_timeseries = zeros(totalNeurons, simulationLength);
u_timeseries = zeros(totalNeurons, simulationLength);
sinusoid_timeseries = zeros(simulationLength,1);
firings = []; % spike timings
timeStep = 1;
timePoint = 1;
for t = 1:timeStep:simulationLength
  
  % Find which neurons have fired (i.e. v>30)
  fired = find(v >= 30);    % indices of spikes
  firings = [firings; t+0*fired,fired]; % Concatenate time and which neuron has fired to 'firings'
  
  % Track v and u
  v_timeseries(:,timePoint) = v;
  u_timeseries(:,timePoint) = u;
  
  % Fulfil reset conditions
  v(fired) = C(fired);
  u(fired) = u(fired) + D(fired);
  
  % Get current synaptic currents (besides the synaptic input, each neuron receives a noisy thalamic input)
  I = vertcat(noise.rs.std .* randn(n_regularSpiking,1) + noise.rs.mean, ...
            noise.fs.std .* randn(n_fastSpiking,1) + noise.fs.mean, ...
            noise.lts.std .* randn(n_lowThreshSpiking,1) + noise.lts.mean);
  I = I + sum(allSynapses(:,fired),2); % Add synaptic weights from neurons that have fired
  alphaWave = sin(t/15.9); % Approximates alpha rhythm
  I = I + alphaWave;
  sinusoid_timeseries(timePoint) = alphaWave;
  
  % Determine current membrane voltage ('v')
  izhikevichEquationForV = (0.04*v.^2) +(5*v) + (140) -(u) +(I);
  update = timeStep * izhikevichEquationForV;
  v = v + update; % Update 'v' 
  
  % Determine current slow recovery status ('u')
  izhikevichEquationForU = A .* ((B.*v) - u);
  u = u + izhikevichEquationForU;
  
  timePoint = timePoint + 1;

end;

%% Plot results

% plot(sinusoid_timeseries)

% Plot graph showing all spikes during simulation window
figure(1); hold on;
excitatorySpikes = find(firings(:,2)<=n_regularSpiking);
inhibitorySpikes = find(firings(:,2)>n_regularSpiking);
plot(firings(excitatorySpikes,1), firings(excitatorySpikes,2), '.', 'Color', 'r');
plot(firings(inhibitorySpikes,1), firings(inhibitorySpikes,2), '.', 'Color', [0 .2 .7]);









%% Reverse connectivity matrices
% connectivity.rs.rs = .0375; % from Regular spiking to...
% connectivity.rs.fs = .125;
% connectivity.rs.lts = .125;
% connectivity.fs.rs = -.25; % from Fast spiking to...
% connectivity.fs.fs = -.15;
% connectivity.fs.lts = -.1;
% connectivity.lts.rs = -.3; % from Low-threshold spiking to...
% connectivity.lts.fs = -.1;
% connectivity.lts.lts = 0;

% %%
%    %% Time specifications:
%    Fs = 1000;                   % samples per second
%    dt = 1/Fs;                   % seconds per sample
%    StopTime = 1;             % seconds
%    t = (0:dt:StopTime-dt)';     % seconds
%    
%    %% Sine wave:
%    Fc = 10;                     % hertz
%    x = sin(2*pi*Fc*t);
%    
%    % Plot the signal versus time:
%    figure;
%    plot(t,x);
%    xlabel('time (in seconds)');
%    title('Signal versus Time');
%    zoom xon;

