function [bs,nave]=data212bi_all(data,segleng,segshift,epleng,maxfreqbin,para);
% usage: bs=data212bi_all(data,segleng,segshift,epleng,maxfreqbin,para)
% 
% calculates cross-bispectra 
% input: 
% data: ndat times nchan matrix each colum is the time-series in one
%             channel;
% segleng: length of each segment in bins, e.g. segleng=1000;  
% segshift: numer of bins by which neighboring segments are shifted;
%           e.g. segshift=segleng/2 makes overlapping segments
% epleng: length of each epoch
% maxfreqbin: max frequency in bins
% para: optional structure:
%       para.segave=0  -> no averaging across segments 
%       para.segave neq 0 -> averaging across segments (default is 0)% \
%       para.subave =1 subtracts the average across epochs,  

%         
% output : 
% 
% bs(m,n,p,f)=<z_m(f)z_n(f)conj(z_p(2f))>
% where F=(f-1)*df is the frequency 

%     
% nave: number of averages


if nargin<6
    para=[];
end




subave=0;
mydetrend=0;
proj=[];
 
   if isfield(para,'detrend')
    mydetrend=para.detrend;
  end 
  
  if isfield(para,'subave')
    subave=para.subave;
  end 
 
  [n nchan]=size(data);
  nep=(n/epleng);
 
 
  if subave==1
      datar=reshape(data,epleng,nep,nchan);
      datam=mean(datar,2);
      for k=1:nep;
          datar(:,k,:)=datar(:,k,:)-datam; 
      end
      data=reshape(datar,n,nchan);
  end
  
 maxfreqbin=floor(min([maxfreqbin,floor(segleng/2)+1])/2);


nseg=floor((epleng-segleng)/segshift)+1; %total number of segments


bs=zeros(nchan,nchan,nchan,maxfreqbin); 


mywindow=repmat(hanning(segleng),1,nchan);
if isfield(para,'mywindow');
    mywindow=repmat(para.mywindow,1,nchan);
end

%figure;plot(mywindow);
nave=0;
for j=1:nep;
    dataep=data((j-1)*epleng+1:j*epleng,:);
    for i=1:nseg; %average over all segments;
        dataloc=dataep((i-1)*segshift+1:(i-1)*segshift+segleng,:);
        if mydetrend==1
           datalocfft=fft(detrend(dataloc,0).*mywindow);
           %datalocfft=datalocfft./abs(datalocfft);
        else
           datalocfft=fft(dataloc.*mywindow);
           %datalocfft=datalocfft./abs(datalocfft);
        end
        
         for f=1:maxfreqbin % for all frequencies
           d1=transpose(datalocfft(f,:));
           d2=conj(datalocfft(2*f-1,:));
           for k=1:nchan;
               bs(:,k,:,f)=bs(:,k,:,f)+reshape(d1*d2*d1(k),nchan,1,nchan);
           end
         end
     nave=nave+1;      
      
    end
   
end

bs=bs/nave;

return;