function [gaussMatrix] = gaussian_random_generator(m, n, class)
gaussMatrix = zeros(m, n, class);
gaussMatrix = random('normal', 0, 1, size(gaussMatrix));
end

