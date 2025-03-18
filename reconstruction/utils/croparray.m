function [ im ] = croparray( img, cropsize )
%CROPARRAY Summary of this function goes here
%   Detailed explanation goes here

if ndims(img) > 2
im = img(1+cropsize(1):end-cropsize(1), 1+cropsize(2):end-cropsize(2), 1+cropsize(3):end-cropsize(3));
else
im = img(1+cropsize(1):end-cropsize(1), 1+cropsize(2):end-cropsize(2));
end

end

