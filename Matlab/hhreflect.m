% Householder reflection
function Y = hhreflect(X) %#codegen

%   Copyright 2011 The MathWorks, Inc.

persistent permutation;

if isempty(permutation)
    permutation = 1:numel(X);
    for i = 1:numel(X)
        ai = int32(rand() * numel(X) + 1);
        bi = int32(rand() * numel(X) + 1);
        if ai > numel(X), ai = int32(numel(X)); end
        if bi > numel(X), bi = int32(numel(X)); end
        [X(ai),X(bi)] = swap(X(ai),X(bi));
    end
end

Y = X(permutation) - 2*sum(X)/numel(X);

function [B,A] = swap(A,B)
