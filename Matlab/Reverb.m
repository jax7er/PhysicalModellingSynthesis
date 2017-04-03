% Implement a reverberation effect using 8 delay lines and householder
% reflection.

%   Copyright 2011 The MathWorks, Inc.

classdef Reverb < handle %#codegen
    % Public parameters for this reverb.
    properties(Access = public)
        Density = 0.5, Feedback = 0.8, Feedthrough = 0.9
    end
    % Private members for this reverb state.
    properties(Access = private)
        NumDelays = 8, SampleRate = 22050, Delay0 = 0.0277,
        delay1, delay2, delay3, delay4, delay5, delay6, delay7, delay8,
        gain1, gain2, gain3, gain4, gain5, gain6, gain7, gain8, out
    end
    
    methods
        % Create a new reverberation effect using given sample rate.
        % Construct the necessary delay lines and gain elements.
        function this = Reverb(SampleRate)
            this.SampleRate = SampleRate; % 44.1 kHz
            d0 = this.Delay0*this.SampleRate;
            d1 = d0*this.Density / this.NumDelays;
            this.delay1 = Delay(get_prime(round(d0 + 1*d1 + d0*rand()*0.1)));
            this.delay2 = Delay(get_prime(round(d0 + 2*d1 + d0*rand()*0.1)));
            this.delay3 = Delay(get_prime(round(d0 + 3*d1 + d0*rand()*0.1)));
            this.delay4 = Delay(get_prime(round(d0 + 4*d1 + d0*rand()*0.1)));
            this.delay5 = Delay(get_prime(round(d0 + 5*d1 + d0*rand()*0.1)));
            this.delay6 = Delay(get_prime(round(d0 + 6*d1 + d0*rand()*0.1)));
            this.delay7 = Delay(get_prime(round(d0 + 7*d1 + d0*rand()*0.1)));
            this.delay8 = Delay(get_prime(round(d0 + 8*d1 + d0*rand()*0.1)));
            this.gain1 = (1/this.NumDelays) + (1/this.NumDelays)*rand()*0.2;
            this.gain2 = (1/this.NumDelays) + (1/this.NumDelays)*rand()*0.2;
            this.gain3 = (1/this.NumDelays) + (1/this.NumDelays)*rand()*0.2;
            this.gain4 = (1/this.NumDelays) + (1/this.NumDelays)*rand()*0.2;
            this.gain5 = (1/this.NumDelays) + (1/this.NumDelays)*rand()*0.2;
            this.gain6 = (1/this.NumDelays) + (1/this.NumDelays)*rand()*0.2;
            this.gain7 = (1/this.NumDelays) + (1/this.NumDelays)*rand()*0.2;
            this.gain8 = (1/this.NumDelays) + (1/this.NumDelays)*rand()*0.2;
        end
                
        % Update internal state of reverb with next input sample
        function update(this, input)
            v = [this.delay1.output() this.delay2.output() this.delay3.output() ...
                 this.delay4.output() this.delay5.output() this.delay6.output() ...
                 this.delay7.output() this.delay8.output()];
            this.out = sum(v) + input*this.Feedthrough;
            g = [this.gain1 this.gain2 this.gain3 this.gain4 ...
                 this.gain5 this.gain6 this.gain7 this.gain8];
            y = (hhreflect(v) .* (1-g)) * this.Feedback + g*input;
            this.delay1.update(y(1));
            this.delay2.update(y(2));
            this.delay3.update(y(3));
            this.delay4.update(y(4));
            this.delay5.update(y(5));
            this.delay6.update(y(6));
            this.delay7.update(y(7));
            this.delay8.update(y(8));
        end
        
        % Get next output sample for this reverb.
        function y = output(this)
            y = this.out;
        end
    end
end
