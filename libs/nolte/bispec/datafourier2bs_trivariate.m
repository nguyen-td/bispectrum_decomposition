function [bs,bsnorm]=datafourier2bs_trivariate(datafft,freqpairs)
% calculates cross-bispectral tensors  from data in fourier domain
%
%
% usage: [bs,bsnorm,nave]=datafourier2bs_fourier(datafft,freqpairs);
%
% input:
% datafft: data in Fourier domain with 3 or 4 indices. First index is
%          frequency, second is channel, 3rd (and 4th) are segement
%          (and epoch). Here, third and fourth index are concatenated and
%          sums are over all segments and epochs.
% freqpairs: pairs of frequency indices.
%            The meaning depends on the physical length of the segments. E.g., if the sampling
%            rate is 1000 Hz and segleng=500, then the segment is 500 ms
%            long, and the frequency resolution is .5 Hz. Then, e.g.
%            f1=6 corresponds to 10 Hz, because counting starts from zero
%            and 10 Hz is the 6.th frequency for that frequency resolution.
%            freqpairs=[6,6] would correspond to the coupling between 10
%            and 20 Hz, i.e. in physical units: f1=10Hz, f2=10Hz and f3=f1+f2=20Hz.
%
%
% output:
% bs: nchan by nchan by nchan by tensor such that
%  bs(i,j,k)=<x(f1)_i*x(f2)_j*conj(x(f1+f2-1)_k)>
%  where f1=freqpairs(f,1) and  f2=freqpairs(f,1),
%  where is the input data in the frequency domain and the average is over
%  epeochs and segments. Note, that the third frequency is set to f1+f2-1,
%  because the first frequency index corresponds to zero, and with
%  f3=f1+f2-1, the physical frequency corresponding to f3 is the sum of the
%  other two frequencies.
%
% bsnnorm: norm to calculate bicoherence using the 3-norm.
%


datafft=datafft(:,:,:);
[nf,nchan,nseg]=size(datafft);

bs=zeros(nchan,nchan,nchan);
bsnorm=zeros(nchan,nchan,nchan);

no=power(mean(abs(datafft).^3,3),1/3);
f1=freqpairs(1);
f2=freqpairs(2);
f3=f1+f2-1;
bsloc=zeros(nchan,nchan,nchan);
for i=1:nseg;
    xx=transpose(datafft(f1,:,i))*datafft(f2,:,i);
    for k=1:nchan;
        bsloc(:,:,k)=xx*conj(datafft(f3,k,i));
    end
    bs=bs+bsloc/nseg;
end


for i=1:nchan;
    for j=1:nchan;
        for k=1:nchan;
            bsnorm(i,j,k)=no(f1,i)*no(f2,j)*no(f3,k);
        end;
    end;
end

    

return;