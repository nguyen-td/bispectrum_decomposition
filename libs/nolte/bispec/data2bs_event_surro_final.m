function [bsall,bsori,nave]=data2bs_event_surro_final(data,segleng,segshift,epleng,freqpairs,para)
% calculates bispectral-tensors  from data in general for event-related
% measurements for original data and surrogates using cycling random shifts
% of epochs for the highest frequency. 
% Based on the data2bs_event function of the METH toolbox by Guido Nolte.
%
% usage: [bsall,bsori,nave]=data2bs_event(data,segleng,segshift,epleng,freqpairs,para);
%
% input: 
% data: ndat times nchan matrix each colum is the time-series in one
%             channel;
% segleng: length of each segment in bins, e.g. segleng=1000;  
% segshift: numer of bins by which neighboring segments are shifted;
%           e.g. segshift=segleng/2 makes overlapping segments
% epleng: leng of each epoch
% freqpairs: pairs of  frequency in bins
% para: otional structure 
% para.nrun   number of randomized data set 
%             default: para.nrun=100
%
% output: 
% bsall: nchan by nchan by nchan by nfreqpairs by number of surrogate data sets
%     tensor such that 
%  bsall(i,j,k,m)=<x(f1)_i*x(f2)_j*conj(x(f1+f2-1)_k)> for the m.th data
%  set, where f1=freqpairs(f,1) and  f2=freqpairs(f,1). 
%  x(f1+f2-1)_k was taken from a random epoch. 
%
%  Note that the first shuffle is NOT the true bispectrum. bsall contains
%  only surrogate bispectra. bs_orig gives the original bispectrum.
%
% bsori: original bispectrum without randomization

% nave: number of averages

[ndat, nchan] = size(data);
[nf, ndum] = size(freqpairs);

% maxfreqbin=sum(freqpairs)-1;
maxfreqbin = sum([max(freqpairs(:, 1)), max(freqpairs(:, 2))]) - 1;

mywindow=repmat(hanning(segleng),1,nchan);
nrun=100; % default number of shuffles
if nargin>5
    if isfield(para,'mywindow');
       mywindow=repmat(para.mywindow,1,nchan);
    end
    if isfield(para,'nrun');
       nrun=para.nrun;
    end
end

if nrun > 0
bsall = zeros(nchan, nchan, nchan, nf, nrun); 
else
    bsall = [];
end





nep=floor(ndat/epleng);
 if nep<10;
        warning('too few epochs for randomization') 
 end
    
nseg=floor((epleng-segleng)/segshift)+1; %total number of segments
datafft=zeros(maxfreqbin,nchan,nseg,nep);

 
%figure;plot(mywindow);

for j=1:nep;
     %disp(j)
    dataep=data((j-1)*epleng+1:j*epleng,:);
    for i=1:nseg; %average over all segments;
      dataloc=dataep((i-1)*segshift+1:(i-1)*segshift+segleng,:);
      datalocfft=fft(detrend(dataloc).*mywindow);
      datafft(:,:,i,j)=datalocfft(1:maxfreqbin,:);
    end
end

bsori = zeros(nchan, nchan, nchan, nf);
fprintf('Progress of %d:', nrun);
for kk=1:nrun+1; 
    if mod(kk, 10) == 0
        fprintf('%d', kk);
    elseif mod(kk, 2) == 0
        fprintf('.');
    end

    nave=0;
    cs=zeros(nchan,nchan,nchan,nf);
    csloc=zeros(nchan,nchan,nchan);
    jcut=ceil((nep-1)*rand)+1;
    jrand=[jcut:nep, 1:jcut-1];
    for j=1:nep
        for i=1:nseg 
            for f=1:nf
                f1 = freqpairs(f,1); f2 = freqpairs(f,2);
                xx=transpose(datafft(f1,:,i,j))*datafft(f2,:,i,j);
                for k=1:nchan
                    if kk==1
                      csloc(:,:,k)=xx*conj(datafft(f1+f2-1,k,i,j));
                    else 
                      csloc(:,:,k)=xx*conj(datafft(f1+f2-1,k,i,jrand(j)));
                    end
                end    
                cs(:,:,:,f) = cs(:,:,:,f) + csloc;
            end
            nave=nave+1;
        end
    end

  cs=cs/nave; 
  if kk==1
       bsori(:,:,:,:)=cs;
  else
       bsall(:,:,:,:,kk)=cs; % nchan, nchan, nchan, nfreq, ishuf
  end
end
bsall = squeeze(bsall(:,:,:,:,2:end));
fprintf('\n');
  

  
    

return;