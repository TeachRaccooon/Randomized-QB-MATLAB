% Symmetric positive definite matrix generator
function [A, S] = gen_spd(n)

Q = randn(n,n);

eigen_mean = 2; 
% can be made anything, even zero 
% used to shift the mode of the distribution
S = diag(abs(eigen_mean+randn(n,1)));
A = Q * S * Q';
return 