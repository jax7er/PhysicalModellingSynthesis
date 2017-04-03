% Implement a delay line using a ring buffer.

%   Copyright 2011 The MathWorks, Inc.

classdef Delay < handle %#codegen
    properties
        % 'buffer' to hold samples (to be delayed)
        % 'index' for the current position in the buffer.
        % head of delay is 'index', tail of delay is 'index+1' modulus
        % buffer size.
        buffer, index;
    end
    methods
        % Construct a new delay line using 'n' number of samples
        % for delay.
        function this = Delay(n)
            this.buffer = zeros(1,n);
            this.index = 1;
        end
        % Update input of the delay line. Recompute internal state.
        function update(this, input)
            this.buffer(this.index) = input;
            this.index = this.index + 1;
            if this.index > numel(this.buffer)
                this.index = 1;
            end
        end
        % Get next output sample.
        function y = output(this)
            y = this.buffer(this.index);
        end
    end
end

