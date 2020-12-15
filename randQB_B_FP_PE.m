
%{
randQB_B_FP_PE - Blocked randomized algorithm for finding QB factorization of a
given matrix A. Iteratively assembles matrices Q and B with a step of
block_size, terminates iterations if either the target rank "k" has been
reached, or the difference berween A and QB, measured with Frobenius norm, 
has gotten smaller than parameter "epsillon". Uses reorthogonalization of 
matrix Q at each iteration to improve accuracy.

Uses Gaussian Random matrix.

This scheme is suitable for matrix A with slowly decaying singular values;
power iterations allow to "strech" out the singular values and achieve 
the result easier.

Is "pass efficient" - esigned for cases with large matrix A so that the amount of passes through
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

    'error' - float - error of approximation A to QB. Here, an error may be
    computed using the following property: |A|^2-|B|^2 = |A - (Q*B)|^2,
    where || denoted a Frobenius norm. 

    'precise_rank' - int - low intrinsic rank of matrix A.

This function is intended to represent a fixed-precision version of an
algorithm, though, an option of specifying target rank is still available.
If the estimate for target rank is unavailable, set parameters 'k' and 's'
such that (k + s) = min(m, n), where m, n are number of rows and columns of
matrix A, respectively. 

Here, we pre-allocate space for Q and B using (k + s) as a parameter and "cut off"
the unnecessary rows and columns if the desired accuracy of approximation
has been reached prior to reaching the specifried traget rank. This seems
like a less expensive way than building Q and B by appending new rows and
columns at each iteration, though it requires having multiplications of
large matrices, mostly consisting of zero elements.

Ref: https://pdfs.semanticscholar.org/99be/879787de8510c099d4a6b1539162b007e4c5.pdf
%}
function [Q, B, error, precise_rank] = randQB_B_FP_PE(A, block_size, epsillon, k, s, power)

    [m, n] = size(A);
    l = k + s;

    norm_A = norm(A, 'fro')^2;
    threshold = epsillon ^ 2;
    termination_flag = false;
    
    %computing a full random matrix
    Omega = gaussian_random_generator(n, k + s);
    
    %power iteartions
    for  j = 1 : power
        [G, ~] = qr(A * Omega, 0);
        [Omega, ~] = qr(transpose(A) * G, 0);
    end
    
    G = A * Omega;
    H = transpose(A) * G;
    
    %Allocating full Q and B
    Q = zeros(m, l);
    B = zeros(l, n);
    
    curr_idx = 1;

    for i = 1 : (l / block_size)
        Omega_i = Omega(:, curr_idx : curr_idx + block_size - 1);
        
        %at step 1, lost of computations are unnecessary
        if i == 1
            Y_i = G(:, curr_idx : curr_idx + block_size - 1);
            [Q_i, R_i] = qr(Y_i, 0);
            B_i = transpose(R_i) \ transpose(H(:, curr_idx : curr_idx + block_size - 1));

            %Inserting new columns and rows into Q and B
            Q(:, (1 : block_size)) = Q_i;
            B((1 : block_size), :) = B_i;
            
            curr_idx = curr_idx + block_size;
        else
            
            Temp = B * Omega_i;
            
            Y_i = G(:, curr_idx : curr_idx + block_size - 1) - (Q * Temp);
            [Q_i, R_] = qr(Y_i, 0);
            
            %reorthogonalization step
            Orthogonalization_buffer = Q_i - (Q * (transpose(Q) * Q_i));
            [Q_i, R_i] = qr(Orthogonalization_buffer, 0);
            R_i = R_i * R_;

            B_i = transpose(R_i) \ (transpose(H(:, curr_idx : curr_idx + block_size - 1)) - (transpose(Y_i) * Q) * B - (transpose(Temp) * B));

            %Inserting new columns and rows into Q and B
            Q(:, curr_idx : curr_idx + block_size - 1) = Q_i;   
            B(curr_idx : curr_idx + block_size - 1, :) = B_i;
            
            curr_idx = curr_idx + block_size;
        end

        norm_B = norm(B_i, 'fro')^2;
        error = norm_A - norm_B;
        
        %precise rank determination
        if error < threshold   
            for j = 1 : block_size
                norm_A = norm_A - norm(B_i(j,:))^2;
                if norm_A < threshold
                    termination_flag = true;
                    break;
                end
            end
        else
            norm_A = error;
        end
        if termination_flag
            precise_rank = (i - 1) * block_size + j;
            break;
        end
    end
    
    if ~termination_flag
        precise_rank = i * block_size;
    end
    
    %Getting rid of extra columns and rows with zero values
    if(i ~= (l / block_size))
        Q = Q(:, 1 : i * block_size);
        B = B(1 : i * block_size, :);
    end
    
    if i == n
        fprintf('Approximation error = %f. Fail to converge within the specified toletance\n\n', sqrt(E));
    end
end

