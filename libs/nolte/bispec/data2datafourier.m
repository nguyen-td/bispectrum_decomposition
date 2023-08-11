function datafft=data2datafourier(data,segleng,segshift,epleng,nf,para)
% calculated data in the Fourier domain. Data are divided into epochs
% and further into (overlapping or not) segments. 
% The output has 4 indices: frequency, channel, segment, epoch. 
% The structure of epochs and segments is similar to the one used in data2cs_event. 
% The data in each segment are windowed (default is Hanning) and Fourier
% transformed. 

[ndat,nchan]=size(data);

mywindow=repmat(hanning(segleng),1,nchan);
if nargin>5
if isfield(para,'mywindow');
    mywindow=repmat(para.mywindow,1,nchan);
end
end

nep=floor(ndat/epleng);
 
    
nseg=floor((epleng-segleng)/segshift)+1; %total number of segments

datafft=zeros(nf,nchan,nseg,nep);

 
%figure;plot(mywindow);

for j=1:nep;
     %disp(j)
    dataep=data((j-1)*epleng+1:j*epleng,:);
    for i=1:nseg; %average over all segments;
      dataloc=dataep((i-1)*segshift+1:(i-1)*segshift+segleng,:);
      datalocfft=fft(detrend(dataloc).*mywindow);
      datafft(:,:,i,j)=datalocfft(1:nf,:);
    end
end

return;