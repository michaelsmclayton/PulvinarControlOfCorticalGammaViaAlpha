function [area1, area2, oscillator1, oscillator2, allAreaFirings] = runModel(params)

    allAreaFirings = cell(params.numberOfTrials, 1);
    timePoints = linspace(1, params.simulationLength, params.simulationLength);

    %% Simulation loop
    for trial = 1:params.numberOfTrials

        disp(['Simulation trial = ' num2str(trial)])

        % Create NeuralAreas
        area1 = NeuralArea(params.simulationLength);
        area2 = NeuralArea(params.simulationLength);

        % Create Oscillators
        oscillator1 = Oscillator(0, params.alphaAmplitude, params.simulationLength);
        oscillator2 = Oscillator(params.phaseShift, params.alphaAmplitude, params.simulationLength);

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
        area1.firings(area1.firings(:,1)<=params.windowToRemove,:) = [];
        area2.firings(area2.firings(:,1)<=params.windowToRemove,:) = [];
        area1.firings(:,1) = area1.firings(:,1) - params.windowToRemove;
        area2.firings(:,1) = area2.firings(:,1) - params.windowToRemove;

        % Store all firings
        allAreaFirings{trial} = {area1.firings, area2.firings};

    end

end