function [img_rsos] = rsos( img, dim )

img_rsos = squeeze(sqrt(sum(img.*conj(img),dim)));

end