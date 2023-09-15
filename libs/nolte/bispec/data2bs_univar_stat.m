
function [bsall,bsallnr,nave]=data2bs_univar_stat(data,segleng,segshift,epleng,maxfreqbins,para)
% calculates bispectrum  from data in general for event-related measurement
% as univariate measures, i.e. always within each sensor, and their null distributions using a shuffling approach.
% Based on the data2bs_event function of the METH toolbox by Guido Nolte.
%
% usage: [cs,csnr,nave]=data2bs_univar(data,segleng,segshift,epleng,maxfreqbins,para);
%
% input:
% data: ndat times nchan matrix each colum is the time-series in one
%             channel;
% segleng: length of each segment in bins, e.g. segleng=1000;
% segshift: numer of bins by which neighboring segments are shifted;
%           e.g. segshift=segleng/2 makes overlapping segments
% epleng: leng of each epoch
% maxfreqbins: maximum frequency in bins, starting at zeros Hertz.
%             The frequency resolution, df, is given by the physical length of a
%             segment, say T. Then df=1/T. E.g. if T=2 seconds, the maxfreqbins=101
%             means that the maximum physical frequency is 50 Hertz.
% para: structure which is eventually used later
%
% output:
% bsall: nchan  by nf by nf by number of surrogate data sets
%     tensor such that
%  bsall(i,f1,f2,m)=<x(f1)_i*x(f2)_j*conj(x(f1+f2-1)_k)> for the m.th data
%  set. x(f1+f2-1)_k was taken from a random epoch.
% First data set is the true bispectrum.
% cs: nchan  by nf by nf tensor for nf frequencies (i.e. nf=maxfreqbins)
%  cs(i,f1,f2)=<x(f1)_i*x(f2)_i*conj(x(f1+f2-1)_i)>
%  where  x is the Fourier-transform of the data of each segment
%
% csn:  corresponding normalization factor defined by
%       csn(i,f1,f2)=N_i(f1) N_i(f2) N_i(f1+f2-1);
%       where N_p(f) is defined as (<abs(x(f)_p)^3>)^(1/3)
%       Bicoherence can be calculated as cs./csn
%
% nave: number of averages

[ndat,nchan]=size(data);
nf=maxfreqbins;
% nf=maxfreqbins-1; % edited by Jasmin on 27.06.2023
nep=floor(ndat/epleng); % total number of epochs
nseg=floor((epleng-segleng)/segshift)+1; %total number of segments

mywindow=repmat(hanning(segleng),1,nchan);
if nargin>5
    if isfield(para,'nrun')
        nrun = para.nrun;
    end
    if isfield(para,'mywindow')
        mywindow=repmat(para.mywindow,1,nchan);
    end
end

% Initialization
if nep<10
    warning('too few epochs for randomization')
end

if nrun>0
    bsall = zeros(nchan,nf,nf,nrun); % bispectrum stat
else
    bsall = [];
end

bsori = zeros(nchan,nf,nf); % bispectrum

fprintf('Progress of %d:', nrun);
for irun = 1:nrun+1
    if mod(irun, 10) == 0
        fprintf('%d', irun);
    elseif mod(irun, 2) == 0
        fprintf('.');
    end

    nave=0;
    cs=zeros(nchan,nf,nf);
    csloc = zeros(nchan,nf,nf);

    if irun == 1 % true bispectrum
        csnr=zeros(nchan,nf,nf);
        csn=zeros(nchan,2*nf-1);
        for j=1:nep
            %disp(j)
            dataep=data((j-1)*epleng+1:j*epleng,:);
            for i=1:nseg %average over all segments;
                dataloc=dataep((i-1)*segshift+1:(i-1)*segshift+segleng,:);
                datalocfft=fft(dataloc.*mywindow);
                datalocfft=datalocfft(1:2*nf-1,:);
                cslocn=((abs(datalocfft)).^3)';
                for ichan=1:nchan
                    xx=hankel(conj(datalocfft(1:2*nf-1,ichan)));
                    csloc(ichan,:,:)=(datalocfft(1:nf,ichan)*transpose(datalocfft(1:nf,ichan))).*xx(1:nf,1:nf);
                    cs(ichan,:,:)=cs(ichan,:,:)+csloc(ichan,:,:);
                end

                nave=nave+1;
                csn=csn+cslocn;
            end
        end

        cs=cs/nave;
        csn=csn/nave;
        csn=power(csn,1/3);
        for i=1:nchan
            for f1=1:nf
                for f2=1:nf
                    csnr(i,f1,f2)=(csn(i,f1)*csn(i,f2)*csn(i,f1+f2-1));
                end
            end
        end
        bsall(:,:,:,irun) = cs; bsallnr(:,:,:) = csnr;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else % Added JDVDV on 05.07.23 to generate surrogate data for stat analysis
        krand = randperm(nep);
        lrand = randperm(nep);
        for j=1:nep
            %disp(j)
            dataep=data((j-1)*epleng+1:j*epleng,:);
            dataep_k = data((krand(j)-1)*epleng+1:krand(j)*epleng,:);
            dataep_l = data((lrand(j)-1)*epleng+1:lrand(j)*epleng,:);

            for i=1:nseg %average over all segments;
                % i)
                dataloc = dataep((i-1)*segshift+1:(i-1)*segshift+segleng,:);
                datalocfft = fft(dataloc.*mywindow);
                datalocfft = datalocfft(1:2*nf-1,:);

                % ii)
                dataloc_k = dataep_k((i-1)*segshift+1:(i-1)*segshift+segleng,:);
                datalocfft_k = fft(dataloc_k.*mywindow);
                datalocfft_k = datalocfft_k(1:2*nf-1,:);

                % iii)
                dataloc_l = dataep_l((i-1)*segshift+1:(i-1)*segshift+segleng,:);
                datalocfft_l = fft(dataloc_l.*mywindow);
                datalocfft_l = datalocfft_l(1:2*nf-1,:);

                for ichan=1:nchan
                    xx=hankel(conj(datalocfft_l(1:2*nf-1,ichan)));
                    csloc(ichan,:,:)=(datalocfft(1:nf,ichan)*transpose(datalocfft_k(1:nf,ichan))).*xx(1:nf,1:nf);
                    cs(ichan,:,:)=cs(ichan,:,:)+csloc(ichan,:,:);
                end
                nave=nave+1;
            end
        end

        cs=cs/nave;
        bsall(:,:,:,irun) = cs; %bsallnr(:,:,:,irun) = csnr;
    end

end

return;

