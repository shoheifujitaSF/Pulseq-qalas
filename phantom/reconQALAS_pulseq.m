addpath(genpath('/cluster/berkin/berkin/Matlab_Code_New/LIBRARY/'));

% Load data
data_file_path=['/autofs/cluster/berkin/fujita/testscan/archive/' ...
    '2023_09_19_bay4_pulseq_traj/meas_MID00465_FID38558_qalas_1x1x3_v4_coil.dat'];

[p,n,e] = fileparts(data_file_path);
basic_file_path=fullfile(p,n);

twix_obj = mapVBVD2(data_file_path);
data_unsorted = twix_obj{end}.image.unsorted();
[adc_len,ncoil,readouts]=size(data_unsorted);

% Load trajectory
traj = readmatrix(['/autofs/cluster/berkin/fujita/testscan/2023_09_19_bay4_pulseq_traj/QALAS_1x1x3mm_R1_TR58_ETL127_revised.txt']);

N = [192 168 56];
nTR = 58;
nETL = 127;
step_size = 5 * nETL;

% Pre-allocate the 4D array
img4d_data = zeros([N(1) N(2) N(3) 5]);

for Acq = 1:5
    indices = [];
    for segment_start = (1+nETL*(Acq-1)):step_size:readouts
        segment_end = min(segment_start + nETL - 1, readouts);
        indices = [indices, segment_start:segment_end];
    end

    data = data_unsorted(:,:,indices(:));
    kspace = zeros([N ncoil]);

    for index=1:(nETL*nTR)
        ky = traj(index,2);
        kz = traj(index,3);
        kspace(:,ky,kz,:) = data(:,:,index);
    end

    im = fft3c2(kspace);
    im3D = abs(sum(im.*conj(im),ndims(im))).^(1/2);

    img4d_data(:,:,:,Acq) = im3D;
end

img4d_nii = make_nii(img4d_data);
save_nii(img4d_nii, [basic_file_path '.nii']);
