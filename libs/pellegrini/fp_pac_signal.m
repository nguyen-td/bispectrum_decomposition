% Generates ground-truth source-level time series with ground-truth univariate 
% (within-region) or bivariate (across-region) interactions and projects
% them to sensor level.
%
% Inputs:
%   params - parameter list
%           params.case - [1, 2, 3] univariate, bivariate or univariate + bivariate
%           params.iInt - [integer] number of interactions
%           params.iReg - [integer] number of voxels per region
%   D      - atlas structure, in this case of the Desikan-Killiany atlas
%
% Outputs:
%   sig           - signal on sensor level
%   brain_noise   - 1/f source level noise projected to sensor level
%   sensor_noise  - sensor noise (white noise)
%   L_save        - leadfield matrix
%   iroi_phase    - ROI indices for phase coupling
%   iroi_amplt    - ROI indices for amplitude coupling
%   D             - atlas structure, in this case of the Desikan-Killiany atlas, see fp_get_Desikan.m for the full documentation
%   fres          - sampling frequency
%   n_trials      - numer of trials/epochs
%   filt          - filter parameters 
%   signal_source - source signals
%
% Copyright (c) 2023 Franziska Pellegrini and Stefan Haufe
% Documentation by Tien Dung Nguyen

function [sig,brain_noise,sensor_noise, L_save,iroi_phase,iroi_amplt,D, fres, n_trials,filt,signal_sources] = fp_pac_signal(params,D)

% set parameters

%total number of samples 
N = 120000; 

%Sampling frequency
fs = 100;
fres = fs; 
frqs = sfreqs(fres, fs); % freqs in Hz

%number of trails and epoch length
n_trials = 60;

%interacting bands 
low = [5 7]; %in Hz
high = [30 34]; %in Hz
band_inds_low = find(frqs >= low(1) & frqs <= low(2)); % indices of interacting low frequencies
band_inds_high = find(frqs >= high(1) & frqs <= high(2)); % indices of interacting high frequencies

%coupling strength = SNR in interacting frequency band 
coupling_snr = 0.6; 


%% filters for band and highpass

[bband_low, aband_low] = butter(5, low/(fs/2));
[bband_high, aband_high] = butter(5, high/(fs/2));
[bhigh, ahigh] = butter(5, 1/(fs/2), 'high');
  
filt.aband_low = aband_low;
filt.bband_low = bband_low;
filt.aband_high = aband_high;
filt.bband_high = bband_high;
filt.ahigh = ahigh;
filt.bhigh = bhigh;
filt.band_inds_low = band_inds_low;
filt.band_inds_high = band_inds_high;
filt.low = low; 
filt.high = high;
 
%% randomly select seed, and in bivariate case also target 

if params.case==1 % in univariate case
    if isempty(params.iroi_phase) && isempty(params.iroi_amplt)
        iroi_phase = randperm(D.nroi,params.iInt)';
        iroi_amplt = [];
    else
        iroi_phase = params.iroi_phase;
        iroi_amplt = [];
    end
    
elseif params.case==2 % in bivariate case 
    if isempty(params.iroi_phase) && isempty(params.iroi_amplt)
        iroi_phase = randperm(D.nroi,params.iInt)';
        iroi_amplt = randperm(D.nroi,params.iInt)';
    else
        iroi_phase = params.iroi_phase;
        iroi_amplt = params.iroi_amplt;
    end

    %be sure that no region is selected twice
    for ii = 1:params.iInt
        while any(iroi_phase==iroi_amplt(ii))
            iroi_amplt(ii) = randi(D.nroi,1,1);
        end
    end
    
elseif params.case==3 % uni + bivariate case (not published) 
    assert(length(params.iInt)==2,...
        'Indicate number of uni- and bivariate interactions in mixed case.')
    iroi_amplt=[];
    if isempty(params.iroi_phase)
        iroi_phase = randperm(D.nroi,sum(params.iInt))'; %select regions for both uni and bivariate interactions
    else
        iroi_phase = params.iroi_phase';
    end
   
    % first entries of iroi_amplt are copies of iroi_phase for uni interactions
    iroi_amplt(1:params.iInt(1)) = iroi_phase(1:params.iInt(1));
    
    % last entries of iroi_amplt are regions for bivar interactions
    if isempty(params.iroi_amplt)
        bivar_a = randperm(D.nroi,params.iInt(2));
    else
        bivar_a = params.iroi_amplt(end-(params.iInt(2)-1):end);
    end

    %be sure that no region is selected twice
    for ii = 1:params.iInt(2)
        while any(iroi_phase==bivar_a(ii))
            bivar_a(ii) = randi(D.nroi,1,1);
        end
    end
    
    iroi_amplt = [iroi_amplt bivar_a]';
    
end

%% indices of signal and noise 

sig_ind = [];
for ii = 1:params.iReg
    if params.case==1 %univariate 
        sig_ind = [sig_ind; (iroi_phase.*params.iReg)-(ii-1)];
        iroi_amplt = iroi_phase; 
        
    elseif params.case==2 %bivariate 
        sig_ind = [sig_ind; (iroi_phase.*params.iReg)-(ii-1), (iroi_amplt.*params.iReg)-(ii-1)];
        
    elseif params.case==3 %uni + bivariate 
        sig_ind = [sig_ind; (iroi_phase.*params.iReg)-(ii-1), (iroi_amplt.*params.iReg)-(ii-1)];
        
    end
end

noise_ind = setdiff(1:params.iReg*D.nroi,sig_ind(:));

%% generate low- and high-frequency signal 

xl = randn(N,  sum(params.iInt)*params.iReg);
for ii = 1:  sum(params.iInt)*params.iReg
    xl(:,ii) = filtfilt(bband_low, aband_low, xl(:,ii));
end

xh = randn(N,  sum(params.iInt)*params.iReg);
for ii = 1: sum(params.iInt)*params.iReg
    xh(:,ii) = filtfilt(bband_high, aband_high, xh(:,ii));
end

%extract phase from low and high signal 
xlh = hilbert(xl);
xlphase = angle(xlh);

xhh = hilbert(xh);
xhphase = angle(xhh);

%ensure that amplitude of high-frequent signal is modulated by phase of
%slow oscillation 
xh = real((1-cos(xlphase)).*exp(1i*xhphase));


%% normalize low and high freq signal to 1/f shape 

xh = xh./norm(xh,'fro');
xl = xl./norm(xl,'fro');

for ii = 1:size(xl,2)
    xl(:,ii) = fp_pinknorm(xl(:,ii));
    xh(:,ii) = fp_pinknorm(xh(:,ii));
end

%% generate interacting sources 

%concenate seed and target voxel activity
if params.case==1 %univariate case 
    %one region contains univariate pac 
    uni_pac = xh + xl;
    %s1 -> N x nInts*nReg
    s1 = uni_pac./norm(uni_pac,'fro'); 
    
elseif params.case==2 %bivariate case 
    %one region contains low signal, the other the modulated high signal 
    %s1 -> N x nInts*2*nReg
    s1 = cat(2,xl,xh);
    s1 = s1./norm(s1(:),'fro');
    
elseif params.case==3 %uni + bivariate case (not published) 
    
    univar_inds = 1:params.iInt(1)*params.iReg;
    bivar_inds = (params.iInt(1)*params.iReg)+1 : sum(params.iInt)*params.iReg;

    %univariate interactions
    uni_pac = xl(:,univar_inds) + xh(:,univar_inds);
    s1_u = uni_pac./norm(uni_pac,'fro');    
    
    %bivariate interactions 
    s1_b = cat(2,xl(:,bivar_inds),xh(:,bivar_inds));
    s1_b = s1_b./norm(s1_b(:),'fro');
    
    %s1 -> N x (nInts(1)+(nInts(2)*2))*nReg
    s1 = cat(2,s1_u,s1_b); 
end

% add pink background noise
backg = mkpinknoise(N, size(s1,2), 1);
backg = backg ./ norm(backg, 'fro');

%combine signal and background noise 
signal_sources = coupling_snr*s1 + (1-coupling_snr)*backg;

%% non-interacting sources

%activity at all voxels but the seed and target voxels
noise_sources = mkpinknoise(N, params.iReg*D.nroi-size(s1,2), 1);

%% leadfield for forward model

L_save = D.leadfield;
L3 = L_save(:, D.sub_ind_cortex, :); % select only voxels that belong to a region 

% multiply with normal direction to get from three to one dipole dimension 
normals = D.normals(D.sub_ind_cortex,:)'; 
for is = 1:numel(D.sub_ind_cortex)
    L_mix(:,is) = squeeze(L3(:,is,:))*squeeze(normals(:,is));
end

%select signal and noise leadfield columns
if params.case==3
    L_sig = L_mix(:,[sig_ind(univar_inds,1)' reshape(sig_ind(bivar_inds,:),[],1)']);
else
    L_sig = L_mix(:,sig_ind(:));
end
L_noise = L_mix(:,noise_ind);


%% project to sensors and generate white noise 

%signal on sensor level 
sig = L_sig * signal_sources';
sig = sig ./ norm(sig(:), 'fro');

%brain noise on sensor level 
try
    brain_noise = L_noise * noise_sources';
    brain_noise = brain_noise ./ norm(brain_noise(:), 'fro');
catch
    error('Something went wrong with seed or target selection.')
end

%white noise on sensor level (sensor noise)
sensor_noise = randn(size(sig));
sensor_noise = sensor_noise ./ norm(sensor_noise(:), 'fro');