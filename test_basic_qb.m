%{
test_basic_qb - Plots the singular values, determined based of rand_QB,
with different values of power iterations parameter (0, 1, 2).
INPUT: initial matrix A; matrix of singular values of A, S; target rank
parameters k & s; number of calls for the algorithms, numiters.
%}

function [] = test_basic_qb(A, S, k, s, numiters)

full_length = k + s;
X = 1:full_length;

title('Accuracy in rank k approximation');
xlabel('k + s');
ylabel('Spectral norm error');

for i = 1:numiters
    
    hold on

    [~, B] = rand_QB(A, k, s, 0);
    S_hat = svd(B);
    S_hat = S_hat(1 : full_length);
    p1 = plot(X, S_hat, 'r');
    p1.Color(4) = 0.1;
    
    
    [~, B_pow] = rand_QB(A, k, s, 1);
    S_hat_pow = svd(B_pow);
    S_hat_pow = S_hat_pow(1 : full_length);
    p2 = plot(X, S_hat_pow, 'g');
    p2.Color(4) = 0.1;
    
    [~, B_pow_2] = rand_QB(A, k, s, 2);
    S_hat_pow_2 = svd(B_pow_2);
    S_hat_pow_2 = S_hat_pow_2(1 : full_length);
    p2 = plot(X, S_hat_pow_2, 'b');
    p2.Color(4) = 0.1;

end

hold on

S = S(1 : full_length);
plot(X, S, 'k');

legend({'Basic RSVD', 'RSVD with p = 1', 'RSVD with p = 2','Exact SVD'},'Location','northeast')

hold off

end

