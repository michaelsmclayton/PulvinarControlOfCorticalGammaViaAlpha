classdef Oscillator < handle % 'handle' allows properties to be updated
    
   %-----------------------------------------------------------------------
   %% INITIAL CLASS PROPERTIES
   %-----------------------------------------------------------------------
   properties
       frequency = 10 / 1000;
       phaseShift; currentVoltage; timeseries; alphaAmplitude;
   end
   
   %-----------------------------------------------------------------------
   %% CLASS METHODS
   %-----------------------------------------------------------------------
   methods
       
       %-------------------------------------------------------------------
       % Constructor method
       %-------------------------------------------------------------------
       function obj = Oscillator(phaseShift, alphaAmplitude, simulationLength)
           obj.phaseShift = phaseShift;
           obj.alphaAmplitude = alphaAmplitude;
           obj.timeseries = zeros(simulationLength,1);
       end
       
       %-------------------------------------------------------------------
       % Update sinusoid
       %-------------------------------------------------------------------
       function update(obj, t) 
           obj.currentVoltage = obj.alphaAmplitude * cos(2*pi * obj.frequency * (t+obj.phaseShift));
           obj.timeseries(t) = obj.currentVoltage;
       end
       
   end
   
end