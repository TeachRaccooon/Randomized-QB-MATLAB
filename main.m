%{
TODO:
Issue with power iterations with blocked_rand_QB_power.
Issue with Frobenius norm trick.
%}

import img_process.*
import gaussian_random_generator.*
import rand_QB.*
import blocked_rand_QB.*

%{
All parameters currently set for test_900x1400.png
If using different matrix object, change k, s, epsillon, block_size and
power manually.
All parameter descriptions are available inside the functions.
%}

filename = 'test_900x1400.png';
A = img_process(filename);

[m, n] = size(A);
epsillon = 0.001;
k = 1000;
s = 250;
block_size = 25;
power = 4;

%[Q, B, error] = rand_QB(A, k, s);
%[Q, B, error] = blocked_rand_QB(A, block_size, epsillon, k, s);
[Q, B, error] = blocked_rand_QB_power(A, block_size, epsillon, k, s, power);
%[Q, B, error] = blocked_rand_QB_large(A, block_size, epsillon, k, s);
%[Q, B, error] = blocked_rand_QB_large_power(A, block_size, epsillon, k, s, power);

%Safer way to measure approximation error
disp(norm(A - (Q * B), 'fro') / norm(A, 'fro'));


