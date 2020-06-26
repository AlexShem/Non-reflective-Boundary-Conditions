function [u, sigma] = GaussianDistrib_3sigma(x, segment)
if nargin == 1
    segment = [x(1) x(end)];
elseif nargin == 2 && length(segment) ~= 2
    error('Length of Domain of Gaussian Distribution should be ')
end
a = segment(1);
b = segment(2);
mn = .5*(a+b);
% Take 3-sigma domain: (-3, 3)
sigma = (b-a) / 6;

u = normpdf(x, mn, sigma);
end
