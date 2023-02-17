function [bsall,bsori,nave]=data2bs_event_surro_final(data,segleng,segshift,epleng,freqpairs,para)
% calculates bispectral-tensors  from data in general for event-related
% measurements for original data and surrogates using cycling random shifts
% of epochs for the highest frequency. 
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
% bsall: nchan by nchan by nchan by number of surrogate data sets
%     tensor such that 
%  bsall(i,j,k,m)=<x(f1)_i*x(f2)_j*conj(x(f1+f2-1)_k)> for the m.th data
%  set. x(f1+f2-1)_k was taken from a random epoch. 
%  where f1=freqpairs(f,1) and  f2=freqpairs(f,1),
%
% bsori: oriiginal bispectrum with randomization

% nave: number of averages

[ndat,nchan]=size(data);
maxfreqbin=sum(freqpairs)-1;

mywindow=repmat(hanning(segleng),1,nchan);
nrun=100; 
if nargin>5
    if isfield(para,'mywindow');
       mywindow=repmat(para.mywindow,1,nchan);
    end
    if isfield(para,'nrun');
       nrun=para.nrun;
    end
end

if nrun>0
bsall=zeros(nchan,nchan,nchan,nrun); 
else
    bsall=[];
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


for kk=0:nrun; 
    nave=0;
    cs=zeros(nchan,nchan,nchan);
    csloc=zeros(nchan,nchan,nchan);
     jcut=ceil((nep-1)*rand)+1;
    jrand=[jcut:nep,1:jcut-1];
  for j=1:nep;
      for i=1:nseg; 
         f1=freqpairs(1);f2=freqpairs(2);
         xx=transpose(datafft(f1,:,i,j))*datafft(f2,:,i,j);
         for k=1:nchan;
             if kk==0;
                  csloc(:,:,k)=xx*conj(datafft(f1+f2-1,k,i,j));
             else 
                  csloc(:,:,k)=xx*conj(datafft(f1+f2-1,k,i,jrand(j)));
             end
         end    
         cs=cs+csloc;
         nave=nave+1;
      end
     
  end

  cs=cs/nave; 
  if kk==0;
       bsori(:,:,:)=cs;
    
  else
       bsall(:,:,:,kk)=cs;  
  end
end
  

  
    

return;