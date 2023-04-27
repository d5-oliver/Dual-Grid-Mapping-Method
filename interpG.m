function G = interpG(xF,xC)
%INTERPG constructs a mapping matrix using linear interpolation.
%
% INTERPG(x,X) assumes that the vector X is a subset of the vector x. This
% requires that the number of elements in the vector x, N, and the length 
% of the vector X, M, are both integers such that (N-1)/(M-1)+1 is also an 
% integer.

N = length(xF);
M = length(xC);
P = (N-1)/(M-1) + 1;

G = zeros(N,M);

for m = 1:M-1

    idx = (m-1)*(P-1)+1:m*(P-1)+1;

    G(idx,m) = (xC(m+1) - xF(idx)) ./ (xC(m+1) - xC(m));
    G(idx,m+1) = (xF(idx) - xC(m)) ./ (xC(m+1) - xC(m));

end

end
