%{
test_time_compare - here, we compare times that it takes for each algorithm
to find all factors of decomposition and "rebuild" the original matrix by
multiplying the factors together. 
%}

function [] = test_time_compare(A, block_size, k, s, numiters)

%S = S(1 : full_length);
X = zeros(1, numiters);
%plot(X, S, 'k');

title('Algorithms execution time comparison');
xlabel('Expected Rank (k + s)');
ylabel('Time (Sec)');

T_B = zeros(1, numiters);
T_PE = zeros(1, numiters);
T_QB = zeros(1, numiters);
T_NY = zeros(1, numiters);

for i = 1:numiters
    
    tic;
    [Q, B] = rand_QB_B_FR(A, block_size, k, s, 0);
    A_k = Q * B;
    T_B(i) = toc;
    
    tic;
    [Q, B] = rand_QB_B_FR_PE(A, block_size, k, s, 0);
    A_k = Q * B;
    T_PE(i) = toc;
    
    tic;
    [Q, B] = rand_QB(A, k, s, 0);
    A_k = Q * B;
    T_QB(i) = toc;
    
    tic;
    [At, ~] = Nyst_Gen(A, k, s);
    A_k = At;
    T_NY(i) = toc;
    
    %rank increase
    X(i) = k;
    k = k + k;
end

hold on

p1 = plot(X, T_B,'color','r');
p1.Marker = 'x';

hold on

p2 = plot(X, T_PE,'color','g');
p2.Marker = 'x';

hold on

p3 = plot(X, T_QB, 'color','b');
p3.Marker = 'x';

hold on

p4 = plot(X, T_NY, 'color','m');
p4.Marker = 'x';


legend({ 'rand\_QB\_B\_FR', 'rand\_QB\_B\_FR\_PE', 'rand\_QB', 'Gen\_Nyst'},'Location','northeast')

hold off

end