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

    'error' - float - error of approximation A to QB. Here, an error is
    computed using the following property: |A|^2-|B|^2 = |A - (Q*B)|^2,
    where || denoted a Frobenius norm. However, it appears that the
    statement holds only if A = QB. In our case, QB only approximates A, so
    it seems that the error is not represented correctly. This statement is
    yet to be revised.

OUTPUT PARAMETERS:
%}
function [Q, B, error] = rand_QB(A, k, s)
    [m, n] = size(A);
    Omega = gaussian_random_generator(n, k + s);
    Q = orth(A * Omega);
    B = transpose(Q) * A;
    
    norm_A = norm(A, 'fro');
    norm_B = norm(B, 'fro');
    norm_A = norm_A * norm_A;
    norm_B = norm_B * norm_B;
    error = norm_A - norm_B;
    
end
