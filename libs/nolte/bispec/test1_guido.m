% addpath ..
% addpath /Volumes/data/matlab/libs/meth/highlevel/

fs = 200;
segleng = 2*fs;
segshift = 2*fs;
epleng = segleng;
maxfreqbins = 200;
frqs=(0:maxfreqbins-1)*fs/segleng;
low = [9 11];
high = [58 62];

nn=1000000;
[xl, xh, pac] = syn_sig(nn, 200, low, high);
[xl2, xh2, pac2] = syn_sig(nn, 200, low, high);
pac=pac+randn(size(pac))*std(pac);
pac2=pac2+randn(size(pac2))*std(pac2);

[cs,csnr,nave] = data2bs_univar(pac, segleng,segshift,epleng,maxfreqbins);

figure; 
subplot(2,2,1);pwelch(pac, 200, [], [], 200);
subplot(2,2,2);imagesc(frqs, frqs, squeeze(abs(cs)))
axis xy
colorbar
subplot(2,2,3);
imagesc(frqs, frqs, squeeze(abs(cs./csnr)))
axis xy
colorbar

% figure; imagesc(squeeze(abs(cs)))
% axis xy


disp('lower sidelobe')


A = eye(2);A(1,2)=2*(rand-.5);A(2,1)=2*(rand-.5);
[cs1,nave] = data2bs_event([pac pac2],segleng,segshift,epleng, [21 101]);

xx1 = cs1 - permute(cs1, [2 1 3]);


[cs2,nave] = data2bs_event([pac pac2]*A,segleng,segshift,epleng, [21 101]);

xx2 = cs2 - permute(cs2, [2 1 3]);


abs([xx2(1,2,1)/cs2(1,1,1),xx2(2,1,1)/cs2(1,1,1),xx2(1,2,2)/cs2(1,1,2),xx2(2,1,2)/cs2(2,2,2)])



%%
[cs2,nave] = data2bs_event([pac pac2]*A,segleng,segshift,epleng, [21 101]);
[cs2norm,nave]=data2bs_threenorm([pac pac2]*A,segleng,segshift,epleng, [21 101]);

bicoh=cs2./cs2norm
bicoha=bicoh-permute(bicoh, [2 1 3])

normsym=(cs2norm+permute(cs2norm,[3 1 2])+permute(cs2norm,[2 3 1]))/3;
b2=cs2./normsym

b2a=b2- permute(b2, [2 1 3])





