function coeffs = fp_fft_coeffs(data,segleng,segshift,epleng,freqpairs)

[ndat,nchan]=size(data);

mywindow=repmat(hanning(segleng),1,nchan);
nep=floor(ndat/epleng);

nseg=floor((epleng-segleng)/segshift)+1; %total number of segments
assert(nseg==1,'only possible with 1 segment')

for j=1:nep
    %disp(j)
    dataloc=data((j-1)*epleng+1:j*epleng,:);
    datalocfft(:,:,j)=fft(detrend(dataloc).*mywindow);
end

coeffs(1,:,:) = datalocfft(freqpairs(1),:,:); %f1
coeffs(2,:,:) = datalocfft(freqpairs(2)-freqpairs(1)+1,:,:); %f2-f1
coeffs(3,:,:) = datalocfft(freqpairs(2),:,:); %f2
coeffs(4,:,:) = datalocfft(freqpairs(2)+freqpairs(1)-1,:,:); %f1+f2