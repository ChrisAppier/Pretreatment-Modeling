%This function generates n data points in a uniform distribution
%between a and b

function [F] = uni_dist(n,a,b)

F = unifrnd(a,b,1,n);

end