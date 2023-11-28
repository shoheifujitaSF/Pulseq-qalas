set(0,'DefaultFigureWindowStyle','docked')
addpath utils/
addpath(genpath('/cluster/berkin/berkin/Matlab_Code_New/LIBRARY/'))
clear;close all;clc;


img_path = '/path_to_nii/qalas_1x1x1mm.nii';
b1_path = '/path_to_b1/b1map.nii';
save_filename = 'prisma_ve11c';
slice_selection = 29;


%--------------------------------------------------------------------------
%% load qalas data  
%--------------------------------------------------------------------------

esp             = 5.8 * 1e-3;
turbo_factor    = 127;
TR          = 4500e-3;
alpha_deg   = 4;
num_reps    = 5;
echo2use    = 1;
gap_between_readouts    = 900e-3;
time2relax_at_the_end   = 0;

img = niftiread(img_path);


%--------------------------------------------------------------------------
%% Load b1 map and resize to qalas
%--------------------------------------------------------------------------

img_b1_load = niftiread(b1_path);
img_b1 = double(img_b1_load)/800;       % reference fa 800
%img_b1 = imrotate(img_b1,90);           %  to match img orientation
img_b1 = imresize3(img_b1, size(img(:,:,:,1))); %  to match img size

%--------------------------------------------------------------------------
%% Select slice
%--------------------------------------------------------------------------

img = img(:,:, slice_selection,:);
img_b1 = img_b1(:,:,slice_selection);

N = [size(img, 1) size(img, 2) size(img, 3)];

% ROI mask for Dicoms
for slc_select = 1:N(3)
   msk(:,:,slc_select) = imfill(rsos(img(:,:,slc_select,:),4) > 0, 'holes'); % or 50(invivo) / 100(nist) thre
end

%--------------------------------------------------------------------------
%% threshold high and low B1 values: use raw b1 map without polyfit
%--------------------------------------------------------------------------

thre_high   = 1.35;
thre_low    = 0.65;

temp        = img_b1 .* msk;

temp(temp > thre_high)  = thre_high;
temp(temp < thre_low)   = thre_low;

temp        = temp .* msk;
img_b1      = temp .* msk;


%--------------------------------------------------------------------------
%% create masks for each b1 value
%--------------------------------------------------------------------------

num_b1_bins = 50; % default: 50

b1_val      = linspace( min(img_b1(msk==1)), max(img_b1(msk==1)), num_b1_bins );
%b1_val      = 1; % no b1 correction
sum_msk     = sum(msk(:));

if length(b1_val) == 1
    % do not use b1 correction
    msk_b1 = msk;
else
    
    msk_b1 = zeross([N,length(b1_val)]);
    
    for t = 1:length(b1_val)
        if t > 1
            msk_b1(:,:,:,t) = (img_b1 <= b1_val(t)) .* (img_b1 > b1_val(t-1));
        else
            msk_b1(:,:,:,t) = msk.*(img_b1 <= b1_val(t));
        end
        
        percent_bin = sum(sum(sum(msk_b1(:,:,:,t),1),2),3) / sum_msk;
        
        if t == length(b1_val)
            msk_b1(:,:,:,t) = img_b1 > b1_val(t-1);
        end
        
        msk_b1(:,:,:,t) = msk_b1(:,:,:,t) .* msk;
    end
end

%--------------------------------------------------------------------------
%% create look up table 
%--------------------------------------------------------------------------
tic
clc
fprintf('generating dictionaries...\n');

t1_entries  = [5:10:3000, 3100:50:5000];
t2_entries  = [1:2:350, 370:20:1000, 1100:100:3000];

T1_entries  = repmat(t1_entries.', [1,length(t2_entries)]).';
T1_entries  = T1_entries(:);
  
T2_entries  = repmat(t2_entries.', [1,length(t1_entries)]);
T2_entries  = T2_entries(:);

t1t2_lut    = cat(2, T1_entries, T2_entries);

% remove cases where T2>T1
idx = 0;
for t = 1:length(t1t2_lut)
    if t1t2_lut(t,1) < t1t2_lut(t,2)
        idx = idx+1;
    end
end

t1t2_lut_prune = zeross( [length(t1t2_lut) - idx, 2] );

idx = 0;
for t = 1:length(t1t2_lut)
    if t1t2_lut(t,1) >= t1t2_lut(t,2)
        idx = idx+1;
        t1t2_lut_prune(idx,:) = t1t2_lut(t,:);
    end
end

disp(['dictionary entries: ', num2str(length(t1t2_lut_prune))])



inv_eff     = 0.7:0.1:1;

signal = zeross([length(t1t2_lut_prune), 5, length(b1_val), length(inv_eff)]);

length_b1_val   = length(b1_val);
length_inv_eff  = length(inv_eff);

cnt             = 0;
iniTime         = clock;
iniTime0        = iniTime;

% parallel computing
delete(gcp('nocreate'))
c = parcluster('local');
total_cores = c.NumWorkers;
parpool(min(ceil(total_cores*.25), length(b1_val)))

parfor b1 = 1:length_b1_val
% for b1 = 1:length_b1_val
    for ie = 1:length_inv_eff
        cnt = cnt + 1;
        
        [Mz, Mxy] = sim_qalas_pd_b1_eff_T2(TR, alpha_deg, esp, turbo_factor, t1t2_lut_prune(:,1)*1e-3, t1t2_lut_prune(:,2)*1e-3, num_reps, echo2use, gap_between_readouts, time2relax_at_the_end, b1_val(b1), inv_eff(ie));
        
        temp = abs(Mxy(:,:,end).');
        
        for n = 1:size(temp,1)
            temp(n,:) = temp(n,:) / sum(abs(temp(n,:)).^2)^0.5;
        end
        
        signal(:,:,b1,ie) = temp;
        
%         if mod(cnt,floor(length_b1_val*length_inv_eff/100)) < 1e-3
%             fprintf('generating dictionary: %.1f%%\n',cnt/(length_b1_val*length_inv_eff)*100);
%             fprintf('elapsed time: %.1f sec\n\n',etime(clock,iniTime));
%             iniTime = clock;
%         end
    end
end
delete(gcp('nocreate'))
fprintf('total elapsed time: %.1f sec\n\n',etime(clock,iniTime0));


%--------------------------------------------------------------------------
%%  concat dictionary across inv_Eff dimension
%--------------------------------------------------------------------------

length_dict = length(t1t2_lut_prune);
dict        = zeross([length_dict * length(inv_eff), 5, length(b1_val)]);

for t = 1:length(inv_eff)
    dict(1 + (t-1)*length_dict : t*length_dict, :, :) = signal(:,:,:,t);
end


%--------------------------------------------------------------------------
%% dictionary fit -> in each slice, bin voxels based on b1 value
%--------------------------------------------------------------------------

fprintf('fitting...\n');

estimate_pd_map = 1;     % set to 1 to estiamte PD map (makes it slower)

T1_map = zeross(N);
T2_map = zeross(N);

PD_map = zeross(N);     % proton density-> won't be estimated if estimate_pd_map is set to 0
IE_map = zeross(N);     % inversion efficiency

iniTime         = clock;
iniTime0        = iniTime;

if length(b1_val) > 1
    % use parfor
    
    delete(gcp('nocreate'))
    c = parcluster('local');    % build the 'local' cluster object
    total_cores = c.NumWorkers; % 48 cores for marduk
    parpool( min(ceil(total_cores*.25), N(3)))
    
    parfor slc_select = 1:N(3)
%     for slc_select = 1:N(3)
        disp(num2str(slc_select))
        
        tic
        for b = 1:length(b1_val)
            msk_slc = msk_b1(:,:,slc_select,b);
            num_vox = sum(msk_slc(:)~=0);
            
            if num_vox > 0
                img_slc = zeross([5, num_vox]);
                
                for t = 1:5
                    temp = sq(img(:,:,slc_select,t));
                    img_slc(t,:) = temp(msk_slc~=0);
                end
                    
                    res = dict(:,:,b) * img_slc; % dot product
                    
                    % find T1, T2 values
                    [~, max_idx] = max(abs(res), [], 1);
                    
                    max_idx_t1t2 = mod(max_idx, length_dict);
                    max_idx_t1t2(max_idx_t1t2==0) = length_dict;
                    
                    res_map = t1t2_lut_prune(max_idx_t1t2,:);
                    
                    % find inv_eff
                    max_idx_ie = 1 + (max_idx - max_idx_t1t2) / length_dict;
                    
                    ie_to_use = inv_eff(max_idx_ie);
                    
                    if estimate_pd_map
                        [Mz_sim, Mxy_sim] = sim_qalas_pd_b1_eff_T2(TR, alpha_deg, esp, turbo_factor, res_map(:,1)*1e-3, res_map(:,2)*1e-3, num_reps, echo2use, gap_between_readouts, time2relax_at_the_end, b1_val(b), ie_to_use.');
                    end
                
                t1_map = zeross(N(1:2));
                t1_map(msk_slc==1) = res_map(:,1);
                
                t2_map = zeross(N(1:2));
                t2_map(msk_slc==1) = res_map(:,2);
                
                ie_map = zeross(N(1:2));
                ie_map(msk_slc==1) = ie_to_use;
                
                if estimate_pd_map
                    Mxy_sim_use = abs(Mxy_sim(:,:,end));
                    
                    scl = zeross([num_vox,1]);
                    
                    for idx = 1:size(Mxy_sim_use,2)
                        scl(idx) = Mxy_sim_use(:,idx) \ img_slc(:,idx);
                    end
                    
                    pd_map = zeross(N(1:2));
                    pd_map(msk_slc~=0) = scl;
                    PD_map(:,:,slc_select) = PD_map(:,:,slc_select) + pd_map;
                end
                
                T1_map(:,:,slc_select) = T1_map(:,:,slc_select) + t1_map;
                T2_map(:,:,slc_select) = T2_map(:,:,slc_select) + t2_map;
                IE_map(:,:,slc_select) = IE_map(:,:,slc_select) + ie_map;
                
            end
        end
        toc
    end
    
    delete(gcp('nocreate'))
    fprintf('total elapsed time: %.1f sec\n\n',etime(clock,iniTime0));

else
    for slc_select = 1:N(3)
        % don't use parfor with b1 map is not used
%         disp(num2str(slc_select))

        for b = 1:length(b1_val)
            msk_slc = msk_b1(:,:,slc_select,b);
            num_vox = sum(msk_slc(:)~=0);

            if num_vox > 0
                img_slc = zeross([5, num_vox]);

                for t = 1:5
                    temp = sq(img(:,:,slc_select,t));
                    img_slc(t,:) = temp(msk_slc~=0);
                end

%                 tic   
                    res = dict(:,:,b) * img_slc; % dot product    

                    % find T1, T2 values    
                    [~, max_idx] = max(abs(res), [], 1);

                    max_idx_t1t2 = mod(max_idx, length_dict);
                    max_idx_t1t2(max_idx_t1t2==0) = length_dict;

                    res_map = t1t2_lut_prune(max_idx_t1t2,:);

                    % find inv_eff                         
                    max_idx_ie = 1 + (max_idx - max_idx_t1t2) / length_dict;

                    ie_to_use = inv_eff(max_idx_ie);
    
                    if estimate_pd_map
                        [Mz_sim, Mxy_sim] = sim_qalas_pd_b1_eff_T2(TR, alpha_deg, esp, turbo_factor, res_map(:,1)*1e-3, res_map(:,2)*1e-3, num_reps, echo2use, gap_between_readouts, time2relax_at_the_end, b1_val(b), ie_to_use.');
                    end
%                 toc

                t1_map = zeross(N(1:2));
                t1_map(msk_slc==1) = res_map(:,1);

                t2_map = zeross(N(1:2));
                t2_map(msk_slc==1) = res_map(:,2);

                ie_map = zeross(N(1:2));
                ie_map(msk_slc==1) = ie_to_use;

                if estimate_pd_map
                    Mxy_sim_use = abs(Mxy_sim(:,:,end));
    
                    scl = zeross([num_vox,1]);
    
                    for idx = 1:size(Mxy_sim_use,2)
                        scl(idx) = Mxy_sim_use(:,idx) \ img_slc(:,idx);
                    end
    
                    pd_map = zeross(N(1:2));
                    pd_map(msk_slc~=0) = scl;     
                    PD_map(:,:,slc_select) = PD_map(:,:,slc_select) + pd_map;
                end

                T1_map(:,:,slc_select) = T1_map(:,:,slc_select) + t1_map;
                T2_map(:,:,slc_select) = T2_map(:,:,slc_select) + t2_map;
                IE_map(:,:,slc_select) = IE_map(:,:,slc_select) + ie_map;
            end
        end
    end
end


%--------------------------------------------------------------------------
%% show results
%--------------------------------------------------------------------------

figure;imshow(T1_map, [0 2500]);colormap hot
figure;imshow(T2_map, [0 400]);colormap hot
% figure;imshow(PD_map, [0 max(PD_map(:))]);colormap hot
% figure;imshow(IE_map, [0.7 1]);

concat = [img(:,:,1), img(:,:,2),img(:,:,3),img(:,:,4),img(:,:,5)];
jimg2(concat)

save(['./fitted_matfiles/',save_filename,'_slice',num2str(slice_selection),'.mat'], 'img', 'img_b1', 'T1_map','T2_map','PD_map','IE_map');