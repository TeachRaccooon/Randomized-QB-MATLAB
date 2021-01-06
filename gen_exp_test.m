%{
gen_exp_test - generates an mxn matrix, whose i-th singular value is
exp(-i/t). Returns the generated matrix A and a vector of its singular
values, s.
%}
function [A, s] = gen_exp_test(m, n, t)
    Buf = randn(m, m);
    [U, ~] = qr(Buf);
    Buf = randn(n, n);
    [V, ~] = qr(Buf);
    p = min(m, n);
    s = zeros(p);
    
    for i = 1:p
        s(i) = exp(-i / t);
    end
    S= spdiags(s', 0, m, n);
    A = U * S * V;