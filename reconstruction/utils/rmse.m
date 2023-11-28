function [ res ] = rmse( img1, img2 )

res = norm(img1(:)-img2(:)) / norm(img2(:));

end

