% Get next prime number 'y' with the closest value to 'x'.
% Return y = x if 'x' is already a prime number.
function y = get_prime(x) %#codegen

%   Copyright 2011 The MathWorks, Inc.

if mod(x,2) == 0
    x = x + 1;
end

y = 0;
while true
    if isprime(x)
        y = x;
        return
    end
    x = x + 2;
end
end

function b = isprime(q)
for i = 2:q-1
    if mod(q,i) == 0
        b = false;
        return
    end
end
b = true;
end
