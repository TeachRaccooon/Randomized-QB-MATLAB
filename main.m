import img_process.*
import gaussian_random_generator.*
import rademacher_random_generator.*

import rand_QB.*
import rand_QB_SP.*
import rand_QB_B_FP.*
import rand_QB_B_FR.*
import rand_QB_B_FP_PE.*
import rand_QB_B_FR_PE.*

import test_basic_qb.*
import test_compare_alg.*
import test_signle_pass.*

% CHOOSE THE TEST MATRIX
%filename = 'test_1500x1000.png';
%A = img_process(filename);
%A = gallery('cauchy', 600); % binomial, chow, dorr, frank, gcdmat, grcar, hanowa, house
%[A, S] = gen_three_test(2000, 2000, 3);
%[A, S] = gen_exp_test(2000, 2000, 75);
%[A, S] = gen_spd(2000);

% SET THE FUNCTION PARAMETERS
epsillon = 0.00000001;
k = 500;
s = 0;
block_size = 500;
power = 0;
numiters = 20;

% RUN TESTS
%test_signle_pass(A, S, k, s, numiters);
%test_basic_qb(A, S, k, s, numiters);
%test_compare_alg(A, S, block_size, k, s, power, numiters);
test_time_compare(A, block_size, k, s, numiters);

