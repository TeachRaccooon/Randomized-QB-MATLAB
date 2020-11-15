function [gaussMatrix] = gaussian_random_generator(m, n)
gaussMatrix = zeros(m, n);
gaussMatrix = random('normal', 0, 1, size(gaussMatrix));
end

