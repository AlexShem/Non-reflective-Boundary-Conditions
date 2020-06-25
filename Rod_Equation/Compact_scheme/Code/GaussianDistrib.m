function [u, sigma] = GaussianDistrib(a, b, x)
% Take 3-sigma domain: (-3, 3)
mn = .5*(a+b);
sigma = (b-a) / 6;
u = normpdf(x, mn, sigma);
end
