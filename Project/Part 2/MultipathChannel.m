L=1000;
N=1;
% How to use:
%
% h = MultipathChannel(L) - generates a vector h of length Lx1 containing L
% channel coefficients for the L paths
% h = MultipathChannel(L) - generates a matrix h of dimention LxN, where
% each column corresponds to L channel coefficients for L paths

h = randn(L,N) + 1i*randn(L,N);

power_profile = exp(-0.5*[0:L-1])';

power_profile = repmat(power_profile,1,N);

h = abs(h).*power_profile;