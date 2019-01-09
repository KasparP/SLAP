function [M,V] = poisson_dev( lambda )
M = 2*exp(-lambda).* (lambda.^(floor(lambda)+1))./factorial(floor(lambda));
V = lambda - M.^2;
end

