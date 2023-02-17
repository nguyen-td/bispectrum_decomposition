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
  
for igroup = 1:2
  
  for isub = 1:Ns(igroup)
    close all
    
    mkdir([fig_dir 'group' num2str(igroup) '_sub' num2str(isub)])
    mkdir([results_dir 'group' num2str(igroup) '_sub' num2str(isub)])
    
    data=textread(allnames{igroup}{isub});

    % univariate bicoherence to check 
    [bs,bsnr,nave]=data2bs_univar(data,segleng,segshift,epleng,maxfreq);

    bsnshow=abs(bs)./bsnr;
    
    [ma in] = max(vec(squeeze(sum(bsnshow))));

    [f1 f2] = ind2sub([51 51], in);
    
    figure; imagesc(squeeze(sum(bsnshow))); colorbar
    export_fig([fig_dir 'group' num2str(igroup) '_sub' num2str(isub) '/bispec_univariate'], ['-r300'], '-a2')
    
    figure; imagesc(squeeze(sum(bsnshow))); colorbar; caxis([0 20])
    export_fig([fig_dir 'group' num2str(igroup) '_sub' num2str(isub) '/bispec_univariate_samescale'], ['-r300'], '-a2')
    
    save([results_dir 'group' num2str(igroup) '_sub' num2str(isub) '/bispec_univariate'], 'bs' , 'bsnr' , 'nave', 'f1', 'f2');
  end
end

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
