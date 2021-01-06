%{
rand_QB_SP - Basic randomized algorithm for finding QB factorization of a
given matrix A, done in a single pass. Hence, power iterations are not
used.

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
%}
function [Q, B, error] = rand_QB_SP(A, k, s)

    class_A = class(A);
    [m, n] = size(A);
    
    Omega = gaussian_random_generator(n, k + s, class_A);
    Omega_ = gaussian_random_generator(m, k + s, class_A);
    
    Y = A * Omega;
    Y_ = transpose(A) * Omega_;
    
    [Q, ~] = qr(Y, 0);
    [Q_, ~] = qr(Y_, 0);  
    
    B = ((transpose(Omega_)*Q)\(transpose(Y_)*Q_))*transpose(Q_);
    
    error = norm(A - (Q * B), 'fro');
end