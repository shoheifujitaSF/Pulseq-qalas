addpath(genpath('/cluster/berkin/berkin/Matlab_Code_New/LIBRARY/'));
clear; clc;

flag.save_nii = 1;
flag.save_kspace = 1;

%--------------------------------------------------------------------------
%% Read twix data and sort kspace based on seq file
%--------------------------------------------------------------------------


data_file_path=['meas_MID00206_FID02149_QALAS_1x1x1mm.dat'];

[p,n,e] = fileparts(data_file_path);
basic_file_path=fullfile(p,n);

twix_obj = mapVBVD2(data_file_path);
data_unsorted = twix_obj{end}.image.unsorted();
[adc_len,ncoil,readouts]=size(data_unsorted);


% Read params from seq file
pulseq_file_path = [p, '/',regexprep(n, 'meas_MID\d+_FID\d+_', ''), '.seq'];
seq=mr.Sequence();
seq.read(pulseq_file_path);

N = seq.getDefinition('Matrix');
nTR = seq.getDefinition('nTR');
nETL = seq.getDefinition('nETL');
%os_factor = seq.getDefinition('os_factor');
os_factor=1;
traj_y = seq.getDefinition('traj_y');
traj_z = seq.getDefinition('traj_z');
step_size = 5 * nETL;

% Pre-allocate the 4D array
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

if (flag.save_nii)
    nii_file_path = fullfile(p,[n,'.nii']);
    img4d_nii = make_nii(img4d_data);
    save_nii(img4d_nii, nii_file_path);
end

if (flag.save_kspace)
    mat_file_path = fullfile(p,[n,'.mat']);
    save(mat_file_path,"kspace","-v7.3");
end
