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

load([results_dir 'bispec_univariate2'], 'ffs', 'ffs_alpha', 'ffs_alphabeta', 'ffs_t');

para.nrun=100; %number of permutations for surrogate, here: 100
      
clear bsall bsori nave
for igroup = 1:2   
  for isub = 1:Ns(igroup)  
    mkdir([results_dir 'group' num2str(igroup) '_sub' num2str(isub)])
    data=textread(allnames{igroup}{isub});
    
    fps = [ffs_alpha{igroup}(isub) ffs_alpha{igroup}(isub); ffs_alpha{igroup}(isub) ffs_alpha{igroup}(isub)*2-1; ffs_t];
    save([results_dir 'group' num2str(igroup) '_sub' num2str(isub) '/bispec_trivariate1'], 'fps');
  
    for ifp = 1:3%size(fps, 1)
      [igroup isub ifp]
      freqpairs = fps(ifp, :);
      tic
      [bsall, bsori, nave] = data2bs_event_surro(data,segleng,segshift,epleng,freqpairs,para);
      toc
      
      save([results_dir 'group' num2str(igroup) '_sub' num2str(isub) '/bispec_trivariate1_fp' num2str(ifp)], 'bsall', 'bsori', 'nave');
    end
  end
end
