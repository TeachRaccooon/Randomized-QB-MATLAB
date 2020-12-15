%{
img_process - small routine for reading in images, transforming them to
grayscale and setting values in range of [0, 10].
%}
function [A] = img_process(filename)
RGB = imread(filename);
A = double(im2gray(RGB));
A = A ./ 255;

%A = uint8(25 * A) - converts back
end

