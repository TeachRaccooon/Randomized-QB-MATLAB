%{
test_basic_qb - Plots the singular values, determined using a single pass
version of rnad_QB.
INPUT: initial matrix A; matrix of singular values of A, S; target rank
parameters k & s; number of calls for the algorithms, numiters.
%}

function [] = test_signle_pass(A, S, k, s, numiters)

full_length = k + s;
X = 1:full_length;

title('Accuracy in rank k approximation');
xlabel('k + s');
ylabel('Spectral norm error');

for i = 1:numiters
    
    hold on
    
    [~, B_SP] = rand_QB_SP(A, k, s);
    S_SP = svd(B_SP);
    S_SP = S_SP(1 : full_length);
    p1 = plot(X, S_SP, 'b');
    p1.Color(4) = 0.1;
    
end

S = S(1 : full_length);
plot(X, S, 'k');

legend({'Single Pass QB', 'Exact SVD'},'Location','northeast')

hold off

end