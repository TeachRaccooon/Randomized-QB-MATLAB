%{
randQB_B_FP_PE - Blocked randomized algorithm for finding QB factorization of a
given matrix A. Iteratively computes matrices Q and B with a step of
block_size, terminates iterations if either the target rank "k + s" has been
reached. Uses reorthogonalization of matrix Q at each iteration to improve accuracy.

This scheme is suitable for matrix A with slowly decaying singular values;
power iterations allow to "strech" out the singular values and achieve 
the result easier.

Is "pass efficient" - esigned for cases with large matrix A so that the amount of passes through
its data is minimized. Hence, requires immediate computation of full random
matrix Omega (using oversampling parameter s, where s is a small integer).

Uses Gaussian Random Projection matrix.

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

This function is intended to represent a fixed-rank version of an
algorithm.

Here, we pre-allocate space for Q and B using (k + s) as a parameter and "cut off"
the unnecessary rows and columns if the desired accuracy of approximation
has been reached prior to reaching the specifried traget rank. This seems
like a less expensive way than building Q and B by appending new rows and
columns at each iteration, though it requires having multiplications of
large matrices, mostly consisting of zero elements.

Ref: https://pdfs.semanticscholar.org/99be/879787de8510c099d4a6b1539162b007e4c5.pdf
%}
function [Q, B, error, precise_rank] = rand_QB_B_FR_PE(A, block_size, k, s, power)
    
    norm_A = norm(A, 'fro');
    if norm_A == 0
        return
    end
    
    precise_rank = 0;

    class_A = class(A);
    [m, n] = size(A);
    max_dim = k + s;
    block_size_init = block_size;
    
    norm_B = 0;
    
    approximation_error = zeros(0, 0, class_A);
    
    Q = zeros(m, 0, class_A);
    B = zeros(0, n, class_A);
    
    %computing a full random matrix
    Omega = gaussian_random_generator(n, k + s, class_A);
    
    %power iteartions
    for  j = 1 : power
        [G, ~] = qr(A * Omega, 0);
        [Omega, ~] = qr(transpose(A) * G, 0);
    end
    
    G = A * Omega;
    H = transpose(A) * G;
    
    curr_idx = 1;

    for i = 1 : (max_dim / block_size)
        Omega_i = Omega(:, curr_idx : curr_idx + block_size - 1);
        
        %at step 1, lost of computations are unnecessary
        if i == 1
            Y_i = G(:, curr_idx : curr_idx + block_size - 1);
            [Q_i, R_i] = qr(Y_i, 0);
            B_i = transpose(R_i) \ transpose(H(:, curr_idx : curr_idx + block_size - 1));

            Q = [Q, Q_i]; %#ok<AGROW>
            B = [B; B_i]; %#ok<AGROW>
            
            curr_idx = curr_idx + block_size;
            B_rows = curr_idx;
        else
            
            Temp = B * Omega_i;
            
            Y_i = G(:, curr_idx : curr_idx + block_size - 1) - (Q * Temp);
            [Q_i, R_] = qr(Y_i, 0);
            
            %reorthogonalization step
            Orthogonalization_buffer = Q_i - (Q * (transpose(Q) * Q_i));
            [Q_i, R_i] = qr(Orthogonalization_buffer, 0);
            R_i = R_i * R_;

            B_i = transpose(R_i) \ (transpose(H(:, curr_idx : curr_idx + block_size - 1)) - (transpose(Y_i) * Q) * B - (transpose(Temp) * B));

            Q = [Q, Q_i]; %#ok<AGROW>
            B = [B; B_i]; %#ok<AGROW>
            
            curr_idx = curr_idx + block_size;
            B_rows = curr_idx;
        end

        % Computing current error approximation with Frobenius norm, and
        % normalizing by norm of A.
        norm_B = hypot(norm_B, norm(B_i, 'fro'));
        approximation_error(i,1) = sqrt(abs(norm_A - norm_B) * (norm_A + norm_B)) / norm_A;
        
        % If the approximation error of the current step is larger than the
        % previous, undo the most current step and terminate execution.
        if (i > 1) && (approximation_error(i) > approximation_error(i-1))
            Q(:, end - block_size + 1 : end) = [];
            B(end - block_size + 1 : end, :) = [];
            approximation_error(i) = [];
            break
        end
        
        
        % If the approximation error did not decrease fast enough in
        % regards to the previous step, increase the block size by its
        % initial value.
        
        %{
        BE CAREFUL WITH EXCEEDING BOUNDS
        if (i > 1) && (approximation_error(i) > approximation_error(i - 1) / 2)
            block_size = min(block_size + block_size_init, max(max_dim - B_rows, 1));
        end
        %}
        
    end
    
    %TODO: this is only for reference - remove. Return approximation_error
    %instead of error.
    error = approximation_error(end);
end

