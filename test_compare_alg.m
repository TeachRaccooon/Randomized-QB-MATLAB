%{
test_compare_alg - Plots the singular values, determined based of rand_QB,
rand_QB_B_FP, rand_QB_B_FP_PE.
INPUT: initial matrix A; matrix of singular values of A, S; target rank
parameters k & s; number of calls for the algorithms, numiters.
%}

function [] = test_compare_alg(A, S, block_size, k, s, power, numiters)

full_length = k + s;

S = S(1 : full_length);
X = 1:full_length;
plot(X, S, 'k');

title('Accuracy in rank k approximation');
xlabel('k + s');
ylabel('Spectral norm error');

for i = 1:numiters
    
    hold on
    
    [~, B_1, ~] = rand_QB_B_FR(A, block_size, k, s, power);
    S_1 = svd(B_1);
    S_1 = S_1(1 : full_length);
    p1 = plot(X, S_1,'color','r');
    p1.Color(4) = 0.1;
    p1.Marker = 'x';
    
    [~, B_2, ~] = rand_QB_B_FR_PE(A, block_size, k, s, power);
    S_2 = svd(B_2);
    S_2 = S_2(1 : full_length);
    p2 = plot(X, S_2,'color','g');
    p2.Color(4) = 0.1;
    p2.Marker = 'x';
    
    
    [~, B_3] = rand_QB(A, k, s, power);
    S_3 = svd(B_3);
    S_3 = S_3(1 : full_length);
    p3 = plot(X, S_3,'color','b');
    p3.Color(4) = 0.1;
    p3.Marker = 'x';
    
    [At_1, ~] = Nyst_Gen(A, k, s);
    S_4 = svd(At_1);
    S_4 = S_4(1 : full_length);
    p4 = plot(X, S_4,'color','m');
    p4.Color(4) = 0.1;
    p4.Marker = 'x';
    
    
end

legend({'Exact SVD', 'rand\_QB\_B\_FR', 'rand\_QB\_B\_FR\_PE', 'rand\_QB', 'Nyst\_Gen'},'Location','northeast')

hold off

end