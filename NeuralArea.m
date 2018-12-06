classdef NeuralArea < handle % 'handle' allows properties to be updated
    
   %% Initial class properties
   properties
       
       % Numbers of each neuron type
       n_regularSpiking = 400;
       n_fastSpiking = 75;
       n_lowThreshSpiking = 25;
       
       % Neuron parameters
       params = struct(...
           'regularSpiking', struct('a', .02, 'b', .2, 'c', -65, 'd', 8), ...
           'fastSpiking', struct('a', .1, 'b', .2, 'c', -65, 'd', 2), ...
           'lowThreshSpiking', struct('a', .02, 'b', .25, 'c', -65, 'd', 2));
       
       % Connectivity parameters
       connectivity = struct(...
           'rs', struct('rs', .0375, 'fs', -.25, 'lts', -.3), ...
           'fs', struct('rs', .125, 'fs', -.15, 'lts', -.1), ...
           'lts', struct('rs', .125, 'fs', -.1, 'lts', 0));
       
       % Noise parameters
       noise = struct(...
           'rs', struct('mean', 4, 'std', 6), ...
           'fs', struct('mean', 5, 'std', 4), ...
           'lts', struct('mean', 4, 'std', 4));
       
       % Parameters to be filled in constructor function
       A; B; C; D; synapses; allSynapses; v; u; totalNeurons; v_timeseries; u_timeseries;
       
       % Contents to be filled during simulation
       fired; firings; I;
 
   end
   
   %% Class methods
   methods
       
       %-------------------------------------------------------------------
       % Constructor method
       %-------------------------------------------------------------------
       function obj = NeuralArea(simulationLength)
           
           % Calculate the total number of neurons
           obj.totalNeurons = obj.n_regularSpiking + obj.n_fastSpiking + obj.n_lowThreshSpiking;
           
           
           % Neuron parameters
           obj.A = vertcat(obj.params.regularSpiking.a * ones(obj.n_regularSpiking,1), ...
                       obj.params.fastSpiking.a * ones(obj.n_fastSpiking,1), ...
                       obj.params.lowThreshSpiking.a * ones(obj.n_lowThreshSpiking,1));
           obj.B = vertcat(obj.params.regularSpiking.b * ones(obj.n_regularSpiking,1), ...
                       obj.params.fastSpiking.b * ones(obj.n_fastSpiking,1), ...
                       obj.params.lowThreshSpiking.b * ones(obj.n_lowThreshSpiking,1));
           obj.C = vertcat(obj.params.regularSpiking.c * ones(obj.n_regularSpiking,1), ...
                       obj.params.fastSpiking.c * ones(obj.n_fastSpiking,1), ...
                       obj.params.lowThreshSpiking.c * ones(obj.n_lowThreshSpiking,1));
           obj.D = vertcat(obj.params.regularSpiking.d * ones(obj.n_regularSpiking,1), ...
                       obj.params.fastSpiking.d * ones(obj.n_fastSpiking,1), ...
                       obj.params.lowThreshSpiking.d * ones(obj.n_lowThreshSpiking,1));
                   
           % Create synaptic connections
           obj.synapses = struct(...
               'rs', struct('rs', obj.connectivity.rs.rs * ones(obj.n_regularSpiking, obj.n_regularSpiking), ...
                            'fs', obj.connectivity.rs.fs * ones(obj.n_regularSpiking, obj.n_fastSpiking), ...
                            'lts', obj.connectivity.rs.lts * ones(obj.n_regularSpiking, obj.n_lowThreshSpiking)), ...
               'fs', struct('rs', obj.connectivity.fs.rs * ones(obj.n_fastSpiking, obj.n_regularSpiking), ...
                            'fs', obj.connectivity.fs.fs * ones(obj.n_fastSpiking, obj.n_fastSpiking), ...
                            'lts', obj.connectivity.fs.lts * ones(obj.n_fastSpiking, obj.n_lowThreshSpiking)), ...
               'lts', struct('rs', obj.connectivity.lts.rs * ones(obj.n_lowThreshSpiking, obj.n_regularSpiking), ...
                            'fs', obj.connectivity.lts.fs * ones(obj.n_lowThreshSpiking, obj.n_fastSpiking), ...
                            'lts', obj.connectivity.lts.lts * ones(obj.n_lowThreshSpiking, obj.n_lowThreshSpiking)));
           obj.allSynapses = vertcat(horzcat(obj.synapses.rs.rs, obj.synapses.rs.fs, obj.synapses.rs.lts), ...
                      horzcat(obj.synapses.fs.rs, obj.synapses.fs.fs, obj.synapses.fs.lts), ...
                      horzcat(obj.synapses.lts.rs, obj.synapses.lts.fs, obj.synapses.lts.lts));
           
           % Set, for all neurons, initial value of V to -65 (membrane potential of the simulated neuron)
           obj.v = -65*ones(obj.totalNeurons,1);
    
           % Set 'u' values as initial 'v' values, multipled by 'b' (sensitivity to 'u')
           obj.u = obj.B.*obj.v;
           
           % Initialise stores for timeseries data
           obj.v_timeseries = zeros(obj.totalNeurons, simulationLength);
           obj.u_timeseries = zeros(obj.totalNeurons, simulationLength);
           
       end
           
       %-------------------------------------------------------------------
       % General update method
       %-------------------------------------------------------------------
       function update(obj, t, timePoint)
           
           % Find which neurons have fired (i.e. v>30)
           obj.fired = find(obj.v >= 30);    % indices of spikes
           obj.firings = [obj.firings; t+0*obj.fired, obj.fired]; % Concatenate time and which neuron has fired to 'firings'
           
           % Record current values of 'v' and 'u'
           obj.v_timeseries(:,timePoint) = obj.v;
           obj.u_timeseries(:,timePoint) = obj.u;
           
           % Fulfil reset conditions
           obj.v(obj.fired) = obj.C(obj.fired);
           obj.u(obj.fired) = obj.u(obj.fired) + obj.D(obj.fired);
            
           % --------------------------------------------------------------
           % Add synaptic currents
           % --------------------------------------------------------------
           
           % Add random noise (with different parameters for each neuron type)
           obj.I = vertcat(obj.noise.rs.std .* randn(obj.n_regularSpiking,1) + obj.noise.rs.mean, ...
                           obj.noise.fs.std .* randn(obj.n_fastSpiking,1) + obj.noise.fs.mean, ...
                           obj.noise.lts.std .* randn(obj.n_lowThreshSpiking,1) + obj.noise.lts.mean);

           % Add inputs from synapses of fired neurons
           obj.I = obj.I + sum(obj.allSynapses(:,obj.fired),2); % Add synaptic weights from neurons that have fired
  
           % Add alpha waves input to inhibitory neurons
           obj.I(obj.n_regularSpiking+1:end) = obj.I(obj.n_regularSpiking+1:end) + 1.5*sin(t/15.9); % Approximates alpha rhythm
           
           % --------------------------------------------------------------
           
           % Determine and update current membrane voltage ('v')
           obj.v = obj.v + 0.5*(0.04*obj.v.^2+5*obj.v+140-obj.u+obj.I); % Two sequential steps of 0.5ms for numerical stability
           obj.v = obj.v + 0.5*(0.04*obj.v.^2+5*obj.v+140-obj.u+obj.I);
  
           % Determine current slow recovery status ('u') 
           obj.u = obj.u + (obj.A .* ((obj.B.*obj.v) - obj.u));
 
       end
   end
end