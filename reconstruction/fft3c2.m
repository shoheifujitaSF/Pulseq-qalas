function res = fft3c2(x)
fctr = size(x,1)*size(x,2)*size(x,3);

res = zeros(size(x));

for m = 1:size(x,5)
    for n=1:size(x,4)
        res(:,:,:,n,m) = 1/sqrt(fctr)*fftshift(fftn(ifftshift(x(:,:,:,n,m))));
    end
end

