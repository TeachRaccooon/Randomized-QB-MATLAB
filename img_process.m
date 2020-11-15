%{
img_process - small routine for reading in images, transforming them to
grayscale and setting values in range of [0, 10].
%}
function [A] = img_process(filename)
RBG = imread(filename);
A = double(im2gray(RBG));
A = A ./ 25;
end

