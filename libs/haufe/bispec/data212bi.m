function [cs1,cs2,dc,dcstd]=data212bi_event(data,segleng,segshift,epleng,maxfreqbin,para);
% usage: [cs1,cs2,dc,dcstd]=data212bi_event(data,segleng,segshift,epleng,maxfreqbin,para)
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
%       para.subave ~= 1 -> no subtraction (default is 1) 
%       IMPORTANT: if you just one epoch (e.g. for continuous data)
%         set para.subave=0 
%         
%       -> averaging across segments (default is 0)
%       para.proj must be a set of vector in channel space,  
%       if it exists then the output raw contains the single trial 
%       Fourier-transform in that channel   
%     
%         
% output for task-related design: 
% cs1: nchan by chan by maxfreqbin  tensor cs1(:,:,f) contains 
%     the 'direct part': 
% cs1(m,n,f)=<x(m,f)x(m,f)conj(x(n,2f))
% 
% cs2 contains the 'indirect part'
% cs2(m,n,f)=<x(m,f)x(n,f)conj(x(m,2f))
%     
% dc=contains the antisymmetric combination 
%    dc(m,n,f)=cs1(m,n,f)-cs2(m,n,f)
%
% dcstd contains the estimated standard deviations of dc. Its real part 
% part is the standard deviation of the real part of dc and its imaginary 
% part is the standard deviation of the imaginary part of dc. 

subave=1; 

if nargin<6
    para=[];
end



segave=0;
mydetrend=0;
proj=[];
  if isfield(para,'segave')
    segave=para.segave;
  end 
   if isfield(para,'detrend')
    mydetrend=para.detrend;
  end 
  if isfield(para,'proj')
    proj=para.proj;
  end 
  if isfield(para,'subave')
    subave=para.subave;
  end 
 
 
 maxfreqbin=floor(min([maxfreqbin,floor(segleng/2)+1])/2);
[ndum,npat]=size(proj);

[ndat,nchan]=size(data);
if npat>0 
   data=data*proj;
   nchan=npat;
end

nep=floor(ndat/epleng);

if subave==1
    dx=reshape(data,epleng,nep,nchan);
    dxm=mean(dx,2);
    for i=1:nep;
        dx(:,i,:)=dx(:,i,:)-dxm;
    end
    data=reshape(dx,epleng*nep,nchan);
end


nseg=floor((epleng-segleng)/segshift)+1; %total number of segments



if segave==0
 cs1=zeros(nchan,nchan,maxfreqbin,nseg); 
 cs2=cs1;
 else
 cs1=zeros(nchan,nchan,maxfreqbin); 
  cs2=cs1;
end

if npat>0
  if segave==0
    cs1=zeros(nchan,nchan,maxfreqbin,nep,nseg); 
    else
    cs1=zeros(nchan,nchan,maxfreqbin,nep);
   end
end
     cs2=cs1;
     dcvarreal=cs1;
      dcvarimag=cs1;


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
   
             if segave==0
                 cs1loc=conj((datalocfft(f,:).^2)'*datalocfft(2*(f-1)+1,:)); 
                 cs1(:,:,f,i)=cs1(:,:,f,i)+cs1loc;
                 cs2loc=transpose(conj((datalocfft(f,:))'*(conj(datalocfft(f,:)).*datalocfft(2*(f-1)+1,:)))); 
                 cs2(:,:,f,i)=cs2(:,:,f,i)+cs2loc;
                  dcvarreal(:,:,f,i)=dcvarreal(:,:,f,i)+(real(cs1loc-cs2loc)).^2;
                 dcvarimag(:,:,f,i)=dcvarimag(:,:,f,i)+(imag(cs1loc-cs2loc)).^2;      
                
            else 
                %disp([i,f,size(datalocfft)])
                 cs1loc=conj((datalocfft(f,:).^2)'*datalocfft(2*(f-1)+1,:)); 
                 cs1(:,:,f)=cs1(:,:,f)+cs1loc; 
                 cs2loc=transpose(conj((datalocfft(f,:))'*(conj(datalocfft(f,:)).*datalocfft(2*(f-1)+1,:)))); 
                 cs2(:,:,f)=cs2(:,:,f)+cs2loc;
                 dcvarreal(:,:,f)=dcvarreal(:,:,f)+(real(cs1loc-cs2loc)).^2;
                 dcvarimag(:,:,f)=dcvarimag(:,:,f)+(imag(cs1loc-cs2loc)).^2;      
           
             end
  
        end
    end
    nave=nave+1;
end

if segave>0;nave=nave*nseg; end;
  cs1=cs1/nave;
  cs2=cs2/nave;
  dc=(cs1-cs2);
  dcvarreal=dcvarreal/(nave);
  dcvarimag=dcvarimag/(nave);
dcstd=sqrt((dcvarreal-(real(dc)).^2)/nave)+sqrt(-1)*sqrt((dcvarimag-(imag(dc)).^2)/nave);


return;