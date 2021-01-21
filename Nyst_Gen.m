%{
Nyst_Gen - plain generalized Nystrom's method - applicable to arbitrary
matrices. 
Based on: https://arxiv.org/pdf/2009.11392.pdf
%}

function [At, error] = Nyst_Gen(A, k, s)

    class_A = class(A);
    [m, n] = size(A);
    
    r = k + s;
    
    X = gaussian_random_generator(n, r, class_A);
    Y = gaussian_random_generator(m, r + (0.5 * r), class_A);
    
    AX = A * X;
    YA = Y' * A;
    [Q, R] = qr(YA * X, 0);

    At = (AX / R) * (Q' * YA);
      
    error = norm(A - At, 'fro');
end

