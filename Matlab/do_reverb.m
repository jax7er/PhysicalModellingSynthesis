function y = do_reverb(y, SampleRate) %#codegen

%   Copyright 2011 The MathWorks, Inc.

% Set definite input types and sizes so we don't need to supply
% example arguments when applying the 'codegen' command.
assert(isa(y, 'double'));
assert(size(y,1) >= 1);
assert(size(y,2) == 1);
assert(isa(SampleRate, 'double'));

persistent reverb;
if isempty(reverb)
    reverb = Reverb(SampleRate);
    reverb.Feedback = 0.95;
end

for i = 1:numel(y)
    reverb.update(y(i));
    y(i) = reverb.output();
end

end
