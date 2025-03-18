function res = ifft3c(x)
fctr = size(x,1)*size(x,2)*size(x,3);
res = zeros(size(x));

for n=1:size(x,4)
    res(:,:,:,n) = sqrt(fctr)*fftshift(ifftn(ifftshift(x(:,:,:,n))));
end
