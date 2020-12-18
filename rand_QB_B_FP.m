%{
rand_QB_B_FP - Blocked randomized algorithm for finding QB factorization of a
given matrix A. Iteratively computes matrices Q and B with a step of
block_size, terminates iterations if either the target rank "k + s" has been
reached, or the difference berween A and QB, measured with Frobenius norm, 
has gotten smaller than parameter "epsillon". Uses reorthogonalization of 
matrix Q at each iteration to improve accuracy.

This scheme is suitable for matrix A with slowly decaying singular values;
power iterations allow to "strech" out the singular values and achieve 
the result easier. Power iteration approach is taken from:
http://people.maths.ox.ac.uk/martinsson/Pubs/2015_randQB.pdf

Uses Gaussian Random Projection matrix.

INPUT PARAMETERS:
    'A' - mat - initial data matrix.

    'block_size' - int - the number of rows of B and columns of Q to be prodced at
    each iteartion.

    'epsillon' - float - tolerance level for the approximation QB to A.

    'k' - int - estimate for the low intrinsic rank of A.

    's' - int - oversampling parameter such that (k + s) <= min(rows(A),
    cols(A)).

    'power' - int - the number of power iterations. If power scheme is not
    used, set parameter to 0.

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

Here, an initial idea was to pre-allocate space for Q and B using (k + s)
as a parameter and "cut off" the unnecessary rows and columns if the desired accuracy 
of approximation has been reached prior to reaching the specifried traget rank. This seems
like a less expensive way than building Q and B by appending new rows and
columns at each iteration, though it requires having multiplications of
large matrices, mostly consisting of zero elements.
An issue arises with dynamic block size, as it becomes unclear as to how
many columns of B and Q to cut off. Hence, currently we are appending data
at each iteration of the loop and suppressing warning meassges.
-TO BE INVESTIGATED
    
Ref: https://pdfs.semanticscholar.org/99be/879787de8510c099d4a6b1539162b007e4c5.pdf
%}

function [Q, B, error, precise_rank] = rand_QB_B_FP(A, block_size, threshold, k, s, power)
    
    precise_rank = 0;

    class_A = class(A);
    [m, n] = size(A);
    max_dim = k + s;
    block_size_init = block_size;
    
    norm_A = norm(A, 'fro');
    norm_B = 0;
    
    if norm_A == 0
        return
    end
    
    approximation_error = zeros(0, 0, class_A);
    
    Q = zeros(m, 0, class_A);
    B = zeros(0, n, class_A);
    
    
    
    % max_dim is never eceeded by dynamic block_size
    for i = 1 : ((max_dim / block_size))
        
        % Computing a small random matrix at each step
        Omega_i = gaussian_random_generator(n, block_size, class_A);
        
        % At step 1, lost of computations are unnecessary
        if i == 1
            %{
            Here and in the similar cases the economy QR is used instead
            of orth(), as for some choices of block_size in relation to
            n, orth() may return smaller amount of columns than block_size
            and the algorithm would fail.
            %}
            [Q_i, ~] = qr(A * Omega_i, 0);
            
            %power iterations
            for j = 1 : power
                [Q_i, ~] = qr(transpose(A) * Q_i, 0);
                [Q_i, ~] = qr(A * Q_i, 0);
            end
            
            B_i = transpose(Q_i) * A;
            
            Q = [Q, Q_i]; %#ok<AGROW>
            B = [B; B_i]; %#ok<AGROW>
            
            B_rows = size(B, 1);
            
        else
            Orthogonalization_buffer = (A * Omega_i) - (Q * (B * Omega_i));
            
            [Q_i, ~] = qr(Orthogonalization_buffer, 0);
            
            %power iterations
            for j = 1 : power
                [Q_i, ~] = qr(transpose(A) * Q_i - transpose(B) * (transpose(Q) * Q_i), 0);
                [Q_i, ~] = qr(A * Q_i - Q * (B * Q_i), 0);
            end
            
            %reorthogonalization step
            Reorthogonalization_buffer = Q_i - (Q * (transpose(Q) * Q_i));
            [Q_i, ~] = qr(Reorthogonalization_buffer, 0);

            B_i = transpose(Q_i) * A - transpose(Q_i) * Q * B;
            
            %Inserting new columns and rows into Q and B
            
            Q = [Q, Q_i]; %#ok<AGROW>
            B = [B; B_i]; %#ok<AGROW>
            
            B_rows = size(B, 1);
            
        end
        
        % Computing current error approximation with Frobenius norm, and
        % normalizing by norm of A.
        norm_B = hypot(norm_B, norm(B_i, 'fro'));
        approximation_error(i,1) = sqrt(abs(norm_A - norm_B) * (norm_A + norm_B)) / norm_A;
        
        % Check for convergence of relative approximation error.
        if approximation_error(i) < threshold       
            break;
        end
        
        
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
        if (i > 1) && (approximation_error(i) > approximation_error(i - 1) / 2)
            block_size = min(block_size + block_size_init, max(max_dim - B_rows, 1));
        end
        
    end
    
    %TODO: this is only for reference - remove.
    error = approximation_error(end);
    
    if i == n
        fprintf('Approximation error = %f. Fail to converge within the specified toletance\n\n', approximation_error(end) / norm_A);
    end
    
end

