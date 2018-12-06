classdef Oscillator < handle % 'handle' allows properties to be updated
    
   %-----------------------------------------------------------------------
   %% INITIAL CLASS PROPERTIES
   %-----------------------------------------------------------------------
   properties
       frequency = 10 / 1000;
       phaseShift; alphaAmplitude; timeseries;
   end
   
   %-----------------------------------------------------------------------
   %% CLASS METHODS
   %-----------------------------------------------------------------------
   methods
       
       %-------------------------------------------------------------------
       % Constructor method
       %-------------------------------------------------------------------
       function obj = Oscillator(phaseShift)
           obj.phaseShift = phaseShift;
       end
       
       %-------------------------------------------------------------------
       % Update sinusoid
       %-------------------------------------------------------------------
       function update(obj, t)   
           obj.alphaAmplitude = 2 * cos(2*pi * obj.frequency * (t+obj.phaseShift));
           obj.timeseries = horzcat(obj.timeseries, obj.alphaAmplitude);
       end
       
   end
   
end