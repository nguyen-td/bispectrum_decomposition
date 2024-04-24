%% Load structs

load('bs_all1.mat')
load('bs_orig1.mat')

%% Antisymmetrization
bs_orig1_anti = bs_orig1 - permute(bs_orig1, [3, 2, 1]); % B_ijk - B_kji
bs_all1_anti = bs_all1 - permute(bs_all1, [3, 2, 1, 4]); % B_ijk - B_kji

%% Compute p-values (over 2nd dimension)
n_shuf = size(bs_all1, 4);
alpha = 0.05;
shuf_vals = squeeze(mean(abs(bs_all1_anti), 2));
% shuf_vals1 = squeeze(mean(abs(bs_all1), 2));
true_val = squeeze(mean(abs(bs_orig1_anti), 2));

P = squeeze(sum(true_val < shuf_vals, 3) ./ n_shuf);
% P = squeeze(sum(true_val < shuf_vals1, 3) ./ n_shuf);
[p_fdr, ~] = fdr(P, alpha);
P_fdr = P;
P_fdr(P > p_fdr) = 1;

P_fdr(P_fdr==0) = 1 / n_shuf;
imagesc(-log10(P_fdr))

% 
imagesc(-log10(squeeze(mean(abs(bs_all1_anti(:,:,:,1)), 2))))
imagesc(-log10(squeeze(shuf_vals(:,:,2))))