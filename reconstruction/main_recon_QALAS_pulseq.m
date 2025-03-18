%--------------------------------------------------------------------------
%% Read twix data and sort kspace based on seq file
%--------------------------------------------------------------------------

addpath(genpath('utils'));

% https://github.com/pulseq/pulseq.git
addpath(genpath('path-to-pulseq/pulseq-master'));

clear; clc;
flag.save_kspace = 0;
flag.save_nii = 0;

load('./data/data_unsorted_nist.mat');

[adc_len,ncoil,readouts]=size(data_unsorted);


% Load params

N = [192 168 56];
nTR = 58;
nETL = 127;
os_factor=1;
step_size = 5 * nETL;

load('./utils/SeqSmplPattern_168x58_R1.mat');



% Pre-allocate
img4d_data = zeros([N(1)*os_factor N(2) N(3) 5]);
kspace = zeros([N(1)*os_factor N(2) N(3) 5 ncoil]);

for contrast = 1:5
    indices = [];
    for segment_start = (1+nETL*(contrast-1)):step_size:readouts
        segment_end = min(segment_start + nETL - 1, readouts);
        indices = [indices, segment_start:segment_end];
    end

    data = data_unsorted(:,:,indices(:));
    kspace_contrast = zeros([N(1)*os_factor N(2) N(3) ncoil]);

    for index=1:(nETL*nTR)
        ky = traj_y(index);
        kz = traj_z(index);
        kspace_contrast(:,ky,kz,:) = data(:,:,index);
    end

    im = fft3c2(kspace_contrast);
    im3D = abs(sum(im.*conj(im),ndims(im))).^(1/2);

    img4d_data(:,:,:,contrast) = im3D;
    kspace(:,:,:,contrast,:) = kspace_contrast;
end



figure(1);
tiledlayout(1,5, 'TileSpacing', 'compact');
for i=1:5
    nexttile;
    imshow(abs(sq(img4d_data(:,:,end/2,i))),[0,4e-4]);title(['Contrast ',num2str(i)])
end

