function res = D(image)

%
% res = D(image)
%
% image = a 2D image
%
% This function computes the finite difference transform of the image
%
% Related functions:
%       adjD , invD 
%
%
% (c) Michael Lustig 2005


Dx = image([2:end,end],:,:) - image;
Dy = image(:,[2:end,end],:) - image;
Dz = image(:,:,[2:end,end]) - image;

res = cat(4,Dx,Dy,Dz);


