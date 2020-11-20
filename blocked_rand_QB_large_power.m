
%{
blocked_rand_QB_large_power - Blocked randomized algorithm for finding QB factorization of a
given matrix A. Iteratively assembles matrices Q and B with a step of
block_size, terminates iterations if either the target rank "k" has been
reached, or the difference berween A and QB, measured with Frobenius norm, 
has gotten smaller than parameter "epsillon". Uses reorthogonalization of 
matrix Q at each iteration to improve accuracy.

Uses Gaussian Random matrix.

Specifically suitable for matrix A with slowly decaying singular values.
Uses power iterations to "strech" them out and achieve the result easier.

Designed for cases with large matrix A so that the amount of passes through
its data is minimized. Hence, requires immediate computation of full random
matrix Omega (using oversampling parameter s, where s is a small integer).

INPUT PARAMETERS:
    'A' - mat - initial data matrix.

    'block_size' - int - the number of rows of B and columns of Q to be prodced at
    each iteartion.

    'epsillon' - float - tolerance level for the approximation QB to A.

    'k' - int - estimate for the low intrinsic rank of A.

    's' - int - oversampling parameter such that (k + s) <= min(rows(A),
    cols(A)).

    'power' - int - the number of power iterations.

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
    yet to be revised. Due to this issue, the '_large' algorithms in this
    repo may not be working correctly, as the current inability to perform
    the "Frobenius norm trick" breaks the idea of minimizing the number of
    passes through the original matrix A. For that reason, it is advised to
    use '_large' algorithms as a fixed-rank solutions rather than a
    fixed-error tolerance.

Here, we pre-allocate space for Q and B using (k + s) as a parameter and "cut off"
the unnecessary rows and columns if the desired accuracy of approximation
has been reached prior to reaching the specifried traget rank. This seems
like a less expensive way than building Q and B by appending new rows and
columns at each iteration, though it requires having multiplications of
large matrices, mostly consisting of zero elements.

Ref: https://pdfs.semanticscholar.org/99be/879787de8510c099d4a6b1539162b007e4c5.pdf
%}
function [Q, B, error] = blocked_rand_QB_large_power(A, block_size, epsillon, k, s, power)

    [m, n] = size(A);
    l = k + s;

    norm_A = norm(A, 'fro');
    norm_A = norm_A * norm_A;
    
    %computing a full random matrix
    Omega = gaussian_random_generator(n, k + s);
    
    %power iteartions
    for  j = 1 : power - 1
        G = orth(A * Omega);
        Omega = orth(transpose(A) * G);
    end
    
    G = A * Omega;
    H = transpose(A) * G;
    
    %Allocating full Q and B
    Q = zeros(m, l);
    B = zeros(l, n);


    for i = 1 : (l / block_size)

        Omega_i = Omega(:, ((i - 1) * block_size + 1 : i * block_size));
        
        %at step 1, lost of computations are unnecessary
        if i == 1
            Y_i = G(:, ((i - 1) * block_size + 1 : i * block_size));
            [Q_i, R_i] = qr(Y_i, 0);
            B_i = transpose(inv(R_i)) * transpose(H(:, ((i - 1) * block_size + 1 : i * block_size)));

            %Inserting new columns and rows into Q and B
            Q(:, (1 : block_size)) = Q_i;
            B((1 : block_size), :) = B_i;
        else
            Y_i = G(:, ((i - 1) * block_size + 1 : i * block_size)) - (Q * (B * Omega_i));
            [Q_i, R_] = qr(Y_i, 0);
            
            %reorthogonalization step
            Orthogonalization_buffer = Q_i - (Q * (transpose(Q) * Q_i));
            [Q_i, R_i] = qr(Orthogonalization_buffer, 0);
            R_i = R_i * R_;

            B_i = transpose(inv(R_i)) * (transpose(H(:, ((i - 1) * block_size + 1 : i * block_size))) - (transpose(Y_i) * Q * B) - ((transpose(Omega_i) * transpose(B)) * B));

            %Inserting new columns and rows into Q and B
            Q(:, ((i - 1) * block_size + 1 : i * block_size)) = Q_i;   
            B(((i - 1) * block_size + 1 : i * block_size), :) = B_i;
        end

        norm_B = norm(B, 'fro');
        norm_B = norm_B * norm_B;
        error = norm_A - norm_B;
        
        if (error < (epsillon * epsillon))
            break
        end
    end
    
    %Getting rid of extra columns and rows with zero values
    if(i ~= (l / block_size))
        Q = Q(:, 1 : i * block_size);
        B = B(1 : i * block_size, :);
    end
end

