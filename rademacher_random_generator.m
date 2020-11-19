function [gaussMatrix] = rademacher_random_generator(m, n)
gaussMatrix = zeros(m, n);
mp = [-1 1]; 
gaussMatrix(:) = mp((rand(1, (m * n))<.5)+1);
end

