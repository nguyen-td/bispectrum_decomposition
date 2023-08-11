addpath(genpath('c:/nolte/meth')); %just to read the data
load data_eeg 
epleng=512;
segleng=256;
segshift=segleng/2;

nf=51; %here, up to 50Hz
 datafft=data2datafourier(data,segleng,segshift,epleng,nf);
 
freqpairs=[11,11]; %nf above must be at least freqpairs(1)+frepairs(2)-1;
                   % here: coupling between 10Hz and 20Hz.  
[bs,bsnorm]=datafourier2bs_trivariate(datafft,freqpairs);
%%
% [nf,nchan,nseg,nep]=size(datafft);
% datafft_bootstrap=datafft(:,:,:,randi(nep,1,nep));
% [bs_bootstrap,bsnorm_bootstrap]=datafourier2bs_trivariate(datafft_bootstrap,freqpairs);
