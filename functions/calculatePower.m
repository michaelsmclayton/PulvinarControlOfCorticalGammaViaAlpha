function [excitatoryPowerStore, inhibitoryPowerStore, fPower, ...
    excitatorySpectrogramStore, inhibitorySpectrogramStore, F] = calculatePower(STM, powerParams)

    addpath(genpath('../../chronux_2_12')) % for Chronux library

    % Initialise stores
    excitatoryPowerStore = [];
    inhibitoryPowerStore = [];
    excitatorySpectrogramStore = [];
    inhibitorySpectrogramStore = [];
    
    % Get frequencies for spectrogram analysis
    F = scal2frq(powerParams.scales,'cmor1-1',1/1000);
    
    % Loop over trials
    for trial = 1:size(STM,2)
        disp(['Analysing power for trial ' num2str(trial)])
        
        % Loop over neural area
        for currentArea = 1:2
            
            % CALCULATE POWER SPECTRA -------------------------------------
            
            % Get excitatory power spectrum
            [excitatoryPowerStore(trial, currentArea, :), fPower, Serr] = ...
                mtspectrumc(STM{trial}.excitatory(currentArea,:), powerParams);
           
            % Get excitatory power spectrum
            [inhibitoryPowerStore(trial, currentArea, :), ~, ~] = ...
                mtspectrumc(STM{trial}.inhibitory(currentArea,:), powerParams);
            
            
            % CALCULATE SPECTROGRAM DATA ----------------------------------
            
            % Get excitatory spectrogram data
            excitatoryCoefficients = cwt(STM{trial}.excitatory(currentArea,:), powerParams.scales,'cmor1-1','plot');
            excitatorySpectrogramStore(currentArea, trial, :, :) = wscalogram('image', squeeze(excitatoryCoefficients), 'scales', F);
            
            % Get inhibitory spectrogram data
            inhibitoryCoefficients = cwt(STM{trial}.inhibitory(currentArea,:), powerParams.scales,'cmor1-1','plot');
            inhibitorySpectrogramStore(currentArea, trial, :, :) = wscalogram('image', squeeze(inhibitoryCoefficients), 'scales', F);
            
            % Normalise spectrogram data
            smoothing = 20;
            timesToAnalyse = 1:1000;
            for f = 1:size(excitatorySpectrogramStore(currentArea, trial, :, :),3)
                % Excitatory data
                excitatorySpectrogramStore(currentArea, trial, f, timesToAnalyse) = excitatorySpectrogramStore(currentArea, trial, f, timesToAnalyse) / max(excitatorySpectrogramStore(currentArea, trial, f, timesToAnalyse));
                excitatorySpectrogramStore(currentArea, trial, f, timesToAnalyse) = smooth(excitatorySpectrogramStore(currentArea, trial, f, timesToAnalyse), smoothing);
                
                % Excitatory data
                inhibitorySpectrogramStore(currentArea, trial, f, timesToAnalyse) = inhibitorySpectrogramStore(currentArea, trial, f, timesToAnalyse) / max(inhibitorySpectrogramStore(currentArea, trial, f, timesToAnalyse));
                inhibitorySpectrogramStore(currentArea, trial, f, timesToAnalyse) = smooth(inhibitorySpectrogramStore(currentArea, trial, f, timesToAnalyse), smoothing);
            end            
 
        end
    end
close all
end

