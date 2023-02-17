% addpath(genpath('d:/nolte/meth'));
addpath /Volumes/data/data/Schizophrenie
% load sa_eeg
% load mri
load allnames

allnames = {allnames_K, allnames_P, allnames_R};

nK=length(allnames_K); %controls
nP=length(allnames_P); %patients
nR=length(allnames_R);  % prodromal  

Ns = [nK nP nR];

fs=256; %sampling rate 
segleng=fs;epleng=2*fs;segshift=segleng/2;
maxfreq=51;

%A=mkfilt_eloreta(sa.V_medium); %eloreta filter for later use
% save A A 
% load A 
%%
%   igroup=1; %healthy group
%   isub=12; %subject 12
    
fig_dir = '../figures/';
results_dir = '../results/';

clear meanbs ffs ffs_alpha ffs_alphabeta
for igroup = 1:2
  
  for isub = 1:Ns(igroup)
    
    load([results_dir 'group' num2str(igroup) '_sub' num2str(isub) '/bispec_univariate'], 'bs' , 'bsnr' , 'nave', 'f1', 'f2');

    bsnshow = abs(bs)./bsnr;
    
    meanbs{igroup}(isub, :, :) = squeeze(mean(bsnshow));
    
    ffs{igroup}(isub, :) = [f1 f2];
    
    alpha_freqs = 8:13;
    [ma in] = max(sum(sum(bsnshow(:, alpha_freqs, alpha_freqs), 3)));
    ffs_alpha{igroup}(isub) = alpha_freqs(in);
    
    beta_freqs = 15:26;
    [ma in] = max(vec(sum(bsnshow(:, alpha_freqs, beta_freqs))));
    [ina inb] = ind2sub([length(alpha_freqs) length(beta_freqs)], in);
    ffs_alphabeta{igroup}(isub, :) = [alpha_freqs(ina) beta_freqs(inb)];

  end
end

[h p ci stats] = ttest2(meanbs{1}, meanbs{2});
figure; imagesc(squeeze(stats.tstat)); colorbar
export_fig([fig_dir 'group' num2str(igroup) '_sub' num2str(isub) '/bispec_univariate_ttest'], ['-r300'], '-a2')
    
% [ma in] = max(vec(stats.tstat));
% [ffst ffst] = ind2sub([51 51], in);

[h p ci stats] = ttest2(ffs_alpha{1}, ffs_alpha{2});

ffs_t = [13 13; 13 22; 2 12]; 

save([results_dir 'bispec_univariate2'], 'ffs', 'ffs_alpha', 'ffs_alphabeta', 'ffs_t');


% clear para;
% para.tunit='Freq/Hz';
% para.funit='Freq/Hz';
% para.timeaxis=[0:50];
% para.freqaxis=[0:50];
% 
% figure;
% showtfinhead(bsnshow,sa.locs_2D,para);
% drawnow
% 
% 
% 
% %%
% 
% alpha_freqs = 8:13;
% beta_freqs = 2*alpha_freqs;
% 
% plot(alpha_freqs, sum(sum(bsnshow(:, alpha_freqs, beta_freqs), 3)));
% 
% [ma in] = max(sum(sum(bsnshow(:, alpha_freqs, beta_freqs), 3)));
% f1 = alpha_freqs(in);
% 
% [ma in] = max(sum(sum(bsnshow, 3)));
% 
% f1=in;
% f2=2*f1-1;  %pick two frequencies (here: 10 Hz and 20 Hz
% 
% %%
% [ma in] = max(vec(squeeze(sum(bsnshow))));
% 
% [f1 f2] = ind2sub([51 51], in);
% 
