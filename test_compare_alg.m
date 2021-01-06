%{
test_compare_alg - Plots the singular values, determined based of rand_QB,
rand_QB_B_FP, rand_QB_B_FP_PE.
INPUT: initial matrix A; matrix of singular values of A, S; target rank
parameters k & s; number of calls for the algorithms, numiters.
%}

function [] = test_compare_alg(A, S, block_size, k, s, numiters)

full_length = k + s;

S = S(1 : full_length);
X = 1:full_length;
plot(X, S, 'k');

title('Accuracy in rank k approximation');
xlabel('k + s');
ylabel('Spectral norm error');

for i = 1:numiters
    
    hold on
    
    [~, B, ~] = rand_QB_B_FR(A, block_size, k, s, 0);
    S_hat = svd(B);
    S_hat = S_hat(1 : full_length);
    p1 = plot(X, S_hat, 'd','color','r','markerfacecolor','r','markersize',8);
    p1.Color(4) = 0.1;
    
    [~, B_pow_1, ~] = rand_QB_B_FR_PE(A, block_size, k, s, 0);
    S_hat_pow_1 = svd(B_pow_1);
    S_hat_pow_1 = S_hat_pow_1(1 : full_length);
    p2 = plot(X, S_hat_pow_1, 'o','color','g','markerfacecolor','b','markersize',9);
    p2.Color(4) = 0.1;
    
    [~, B_pow_2] = rand_QB(A, k, s, 0);
    S_hat_pow_2 = svd(B_pow_2);
    S_hat_pow_2 = S_hat_pow_2(1 : full_length);
    p2 = plot(X, S_hat_pow_2, 'x','color','b','markerfacecolor','b','markersize',9);
    p2.Color(4) = 0.1;
    
end

legend({'Exact SVD', 'rand\_QB\_B\_FR', 'rand\_QB\_B\_FR\_PE', 'rand\_QB'},'Location','northeast')

hold off

end