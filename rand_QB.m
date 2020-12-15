%{
rand_QB - Basic randomized algorithm for finding QB factorization of a
given matrix A. 

INPUT PARAMETERS:
    'A' - mat - initial data matrix.

    'k' - int - estimate for the low intrinsic rank of A.

    's' - int - oversampling parameter such that (k + s) <= min(rows(A),
    cols(A)).

OUTPUT PARAMETERS:
    'Q' - mat - orthonormal matrix of with the number of rows same as A and
    number of columns <= (k + s). A is approximately equal to QB.

    'B' - mat - matrix of with the number of columns same as A and
    number of rows <= (k + s). Serves as a "small approximation" for A. For
    instance, in randomized Singular Value Decomposition, to find an approximation 
    for an SVD of A, the function SVD would be applied to matrix B (Then,
    matrix U shall be multiplied by Q on the left to complete the
    procedure).

    'error' - float - error of approximation A to QB. Here, an error may be
    computed using the following property: |A|^2-|B|^2 = |A - (Q*B)|^2,
    where || denoted a Frobenius norm. 

OUTPUT PARAMETERS:
%}
function [Q, B, error] = rand_QB(A, k, s)
    [~, n] = size(A);
    Omega = gaussian_random_generator(n, k + s);
    [Q, ~] = qr(A * Omega, 0);
    B = transpose(Q) * A;
    
    error = norm(A - (Q * B), 'fro');
end

