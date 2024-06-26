function [ res ] = svd_apply3d( in, comp_mtx )
% svd coil compression for 3d data
% assumes that coil axis is the 4th dimension

mtx_size = size(in(:,:,:,1));
 
temp = reshape(in, prod(mtx_size), []);

res = reshape(temp * comp_mtx, [mtx_size, size(comp_mtx,2)]);

end