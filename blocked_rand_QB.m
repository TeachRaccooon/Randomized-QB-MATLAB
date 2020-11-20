%{
blocked_rand_QB - Blocked randomized algorithm for finding QB factorization of a
given matrix A. Iteratively computes matrices Q and B with a step of
block_size, terminates iterations if either the target rank "k + s" has been
reached, or the difference berween A and QB, measured with Frobenius norm, 
has gotten smaller than parameter "epsillon". Uses reorthogonalization of 
matrix Q at each iteration to improve accuracy.

Uses Gaussian Random matrix.

INPUT PARAMETERS:
    'A' - mat - initial data matrix.

    'block_size' - int - the number of rows of B and columns of Q to be prodced at
    each iteartion.

    'epsillon' - float - tolerance level for the approximation QB to A.

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

    'error' - float - error of approximation A to QB. Currently, measured
    by |A - QB|, where || denoted a Frobenius norm. Here, an error may be
    computed using the following property: |A|^2-|B|^2 = |A - (Q*B)|^2,
    where || denoted a Frobenius norm. However, it appears that the
    statement holds only if A = QB. In our case, QB only approximates A, so
    it seems that the error is not represented correctly. This statement is
    yet to be revised. 

Here, we pre-allocate space for Q and B using (k + s) as a parameter and "cut off"
the unnecessary rows and columns if the desired accuracy of approximation
has been reached prior to reaching the specifried traget rank. This seems
like a less expensive way than building Q and B by appending new rows and
columns at each iteration, though it requires having multiplications of
large matrices, mostly consisting of zero elements.
    

Ref: https://pdfs.semanticscholar.org/99be/879787de8510c099d4a6b1539162b007e4c5.pdf
%}
function [Q, B, error] = blocked_rand_QB(A, block_size, epsillon, k, s)

    [m, n] = size(A);
    l = k + s;
    
    %Allocating full Q and B
    Q = zeros(m, l);
    B = zeros(l, n);

    for i = 1 : ((l / block_size))
        
        %computing a small random matrix at each step
        Omega_i = gaussian_random_generator(n, block_size);
        
        %at step 1, lost of computations are unnecessary
        if i == 1
            %{
            Here and in the similar cases the economy QR is used instead
            of orth(), as for some choices of block_size in relation to
            n, orth() may return smaller amount of columns than block_size
            and the algorithm would fail.
            %}
            [Q_i, ~] = qr(A * Omega_i, 0);
            
            B_i = transpose(Q_i) * A;
            
            %Inserting new columns and rows into Q and B
            Q(:, (1 : block_size)) = Q_i;
            B((1 : block_size), :) = B_i;
        else
            Orthogonalization_buffer = (A * Omega_i) - (Q * (B * Omega_i));
            [Q_i, ~] = qr(Orthogonalization_buffer, 0);
            
            %reorthogonalization step
            Reorthogonalization_buffer = Q_i - (Q * (transpose(Q) * Q_i));
            [Q_i, ~] = qr(Reorthogonalization_buffer, 0);

            B_i = transpose(Q_i) * A;
            
            %Inserting new columns and rows into Q and B
            Q(:, ((i - 1) * block_size + 1 : i * block_size)) = Q_i;   
            B(((i - 1) * block_size + 1 : i * block_size), :) = B_i;
        end

        error = norm(A - (Q * B), 'fro');
        
        if (error < epsillon)
            break
        end
    end

    %Getting rid of extra columns and rows with zero values
    if(i ~= (l / block_size))
        Q = Q(:, 1 : i * block_size);
        B = B(1 : i * block_size, :);
    end
end

