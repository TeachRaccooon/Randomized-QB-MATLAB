%{
gen_three_test - generates one of the three random test matrices, 
based on the version parameter.
%}

function [A, s] = gen_three_test(m, n, version)
    Buf = randn(m, m);
    [U, ~] = qr(Buf);
    Buf = randn(n, n);
    [V, ~] = qr(Buf);
    p = min(m,n);
    s = zeros(1, p);
    
    switch version
        case 1 
            s = (1:p).^(-2);
        case 2 
            s = exp(-(1:p)/7);
        case 3 
            s = 0.0001+1./(1+exp((1:p)-30));
        otherwise
            disp('Invalid number. Use [1,2,3].');
            exit;
    end

    S= spdiags(s', 0, m, n);
    A = U * S * V;