% addpath(genpath('d:/nolte/meth'));
addpath /Volumes/data/data/Schizophrenie
load sa_eeg
load mri
load allnames

nK=length(allnames_K); %controls
nP=length(allnames_P); %patients
nR=length(allnames_R);  % prodromal  

fs=256; %sampling rate 
segleng=fs;epleng=2*fs;segshift=segleng/2;
maxfreq=51;

%A=mkfilt_eloreta(sa.V_medium); %eloreta filter for later use
% save A A 
load A 
%%
    igroup=1; %healthy group
    isub=12; %subject 12
    
fig_dir = '../figures/';
results_dir = '../results/';
  
for igroup = 1:3
  for isub = 1:12
    
    if igroup==1;
       data=textread(allnames_K{isub});
    else
       data=textread(allnames_P{isub});
    end

    % univariate bicoherence to check 
    [bs,bsnr,nave]=data2bs_univar(data,segleng,segshift,epleng,maxfreq);

    bsnshow=abs(bs)./bsnr;
    
    [ma in] = max(vec(squeeze(sum(bsnshow))));

    [f1 f2] = ind2sub([51 51], in);
    
    figure; imagesc(squeeze(sum(bsnshow)));

    export_fig([fig_dir 'bispec_univariate'], ['-r300'], '-a2')
  end
end

keyboard

clear para;
para.tunit='Freq/Hz';
para.funit='Freq/Hz';
para.timeaxis=[0:50];
para.freqaxis=[0:50];

figure;
showtfinhead(bsnshow,sa.locs_2D,para);
drawnow



%%

alpha_freqs = 8:13;
beta_freqs = 2*alpha_freqs;

plot(alpha_freqs, sum(sum(bsnshow(:, alpha_freqs, beta_freqs), 3)));

[ma in] = max(sum(sum(bsnshow(:, alpha_freqs, beta_freqs), 3)));
f1 = alpha_freqs(in);

[ma in] = max(sum(sum(bsnshow, 3)));

f1=in;
f2=2*f1-1;  %pick two frequencies (here: 10 Hz and 20 Hz

%%
[ma in] = max(vec(squeeze(sum(bsnshow))));

[f1 f2] = ind2sub([51 51], in);

%%

freqpairs=[f1,f2];

% calculate full bispectrum at chosen frequencies 
para.nrun=0; %number of permutations for surrogate, here: none
[bsall,bsori,nave]=data2bs_event_surro(data,segleng,segshift,epleng,freqpairs,para);

bsori = bsori - permute(bsori, [1 3 2]) - permute(bsori, [2 1 3]) - permute(bsori, [3 1 2]);


%%
% make a low dimensional fit: 
n=3 % number of sources  
 [a,d,err1,err_all,bsmodel]=bsfit(bsori,n);

 
 for ii = 1:n
   figure; showfield(a(:, ii), sa.locs_2D);
 end
 
 %%

 % Make a scan with SC-MUSIC 
 ns=3;
 [s_all,vmax_all,imax_all,dips_mom_all,dips_loc_all]=sc_music(a,sa.V_medium,ns,sa.grid_medium);
 % 
 
 % View sources 
for k=1:ns
 z=1./(1-s_all(:,k));
 zz=[sa.grid_medium,z];
 figure;showmri_transp(mri,[],zz);
 drawnow
 
end

%%

  % Make an inverse with MOCA 
  % First, find mixed  sources 
  [nchan nvoxel ndum]=size(A);
  F=zeros(nvoxel,ndum,n); % intitialze distributions for n sources 
  
  for i=1:n; 
    for k=1:ndum
      F(:,k,i)=A(:,:,k)'*a(:,i); 
    end;
  end
  
 
  % now unmix the sources 
  [Fout,wall]=moca_ncomp(F);
  
  
 % view the result
for i=1:n
    Floc=sqrt(sum(Fout(:,:,i).^2,2));
    Floc(Floc < prctile(Floc, 80)) = 0;
 zz=[sa.grid_medium,Floc];
 figure;showmri_transp(mri,[],zz);
 drawnow
 
end

%%

a_moca = a*wall';
 
for ii = 1:n
 figure; showfield(a_moca(:, ii), sa.locs_2D);
end
