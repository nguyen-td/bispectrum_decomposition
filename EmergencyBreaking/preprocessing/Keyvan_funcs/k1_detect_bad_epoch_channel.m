
function [reject] = k1_detect_bad_epoch_channel(X,par)

%This function detects noisy channels and epoch

%   X: EEG data in the form of: Channels * Time * Epoch
%      if X is in EEGLAB data format set par = []
%
%   par: if X is not in EEGLAB format, sampling rate must be set:
%   par. srate = 100 Hz % example

%   Optional inputs:






if ~isfield(par, 'deviationAnalysis')
    par.deviationAnalysis.run = 1;
end

if ~isfield(par.deviationAnalysis, 'threshold')
    par.deviationAnalysis.threshold = 4.5;
end


if ~isfield(par, 'HFAnalysis')
    par.HFAnalysis.run = 1;
     par.HFAnalysis.channelRejection = 1;
end
par.HFAnalysis.threshold = 5;   % Noise threshold
par.HFAnalysis.lpFilterType = 'low';
par.HFAnalysis.lpFilterRange = 35;



par.epoch.chanDeviation = [];
par.epoch.chanExDeviation = [];
par.epoch.epochDeviation = [];
par.epoch.epochExDeviation = [];
par.chan.chanDeviation = [];
par.epoch.chanHFNoise = [];
par.epoch.epochHFNoise = [];


par.epoch.chanDeviation2 = [];
par.epoch.chanExDeviation2 = [];
par.epoch.epochDeviation2 = [];
par.epoch.epochExDeviation2 = [];
par.epoch.chanHFNoise2 = [];
par.epoch.epochHFNoise2 = [];

par.chanRejectFinal = [];
par.epochRejectFinal = [];

reject = par;




if isstruct(X)
    srate = X.srate;
    X = X.data;
else
    srate = par.srate;
end

if ndims(X)~=3
    error('Convert data to 3D format: Channels * Time * Epochs! ')
end

[Nchan, Ntp, Nepoch] = size(X);

%% %% BAD CHANNEL DETECTION



if par.deviationAnalysis.run
    
    %% METHOD 1.1:  Deviation across channels; epoch by epoch
    reject.epoch.zValChanDeviation = zeros(Nchan, Nepoch);
    for iepoch = 1:Nepoch
        chanDeviation = 0.7413 * iqr(squeeze(X(:,:,iepoch)), 2); % Robust estimate of SD
        chanDeviationSD =  0.7413 * iqr(chanDeviation);
        chanDeviationMedian = nanmedian(chanDeviation);
        reject.epoch.zValChanDeviation(:,iepoch) = ...
            (chanDeviation - chanDeviationMedian) / chanDeviationSD;
    end
    reject.epoch.chanDeviation = (abs(reject.epoch.zValChanDeviation) > reject.deviationAnalysis.threshold);
    reject.epoch.chanExDeviation = (abs(reject.epoch.zValChanDeviation) > 2*reject.deviationAnalysis.threshold);
    
    
    % TO BE TESTED MORE
    %% METHOD 1.2:  Deviation across epochs; channel by channel
    reject.epoch.zValEpochDeviation = zeros(Nchan, Nepoch);
    for ichan = 1:Nchan
        epochDeviation = 0.7413 * iqr(squeeze(X(ichan,:,:)), 1); % Robust estimate of SD
        epochDeviationSD =  0.7413 * iqr(epochDeviation);
        epochDeviationMedian = nanmedian(epochDeviation);
        reject.epoch.zValEpochDeviation(ichan,:) = ...
            (epochDeviation - epochDeviationMedian) / epochDeviationSD;
    end
    reject.epoch.epochDeviation = (abs(reject.epoch.zValEpochDeviation) > reject.deviationAnalysis.threshold);
    reject.epoch.epochExDeviation = (abs(reject.epoch.zValEpochDeviation) > 2 * reject.deviationAnalysis.threshold);
    
    
    
    %% METHOD 1.3:  Deviation across channels; all time points (X is transformed to 2D form: Nchannel * (Ntp*Nepoch)
    reject.chan.zValChanDeviation = zeros(Nchan);
    chanDeviation = 0.7413 *iqr(reshape(X,Nchan,[]), 2); % Robust estimate of SD
    chanDeviationSD =  0.7413 * iqr(chanDeviation);
    chanDeviationMedian = nanmedian(chanDeviation);
    reject.chan.zValChanDeviation = ...
        (chanDeviation - chanDeviationMedian) / chanDeviationSD;
    reject.chan.chanDeviation = (abs(reject.chan.zValChanDeviation) > reject.deviationAnalysis.threshold);
    
end





if par.HFAnalysis.run
    
    [b, a] = butter(2, par.HFAnalysis.lpFilterRange/(srate/2), par.HFAnalysis.lpFilterType);
    Xlp = reshape(filtfilt(b, a, reshape(double(X), Nchan, [])')',Nchan,Ntp,Nepoch);
    
    
    %% METHOD 2.1:  High Frequency (HF) Noise detection: across channels; epoch by epoch
    reject.epoch.zValChanHFNoise = zeros(Nchan, Nepoch);
    for iepoch = 1:Nepoch
        noisiness = mad(X(:,:,iepoch)- Xlp(:,:,iepoch), 1, 2)./mad(Xlp(:,:,iepoch), 1, 2);
        noisinessMedian = nanmedian(noisiness);
        noisinessSD = mad(noisiness, 1, 1)*1.4826; % for a normal distribution sigma=1.4826*MAD
        reject.epoch.zValChanHFNoise(:,iepoch) = (noisiness - noisinessMedian) ./ noisinessSD;
    end
    reject.epoch.chanHFNoise = reject.epoch.zValChanHFNoise >  reject.HFAnalysis.threshold;
    
    
    
    %% METHOD 2.2:  High Frequency (HF) Noise detection: across epochs; channel by channel
    reject.epoch.zValEpochHFNoise = zeros(Nchan, Nepoch);
    for ichan = 1:Nchan
        noisiness = mad(squeeze(X(ichan,:,:)- Xlp(ichan,:,:)), 1, 1)./mad(squeeze(Xlp(ichan,:,:)), 1, 1);
        noisinessMedian = nanmedian(noisiness);
        noisinessSD = mad(noisiness, 1, 2)*1.4826; % for a normal distribution sigma=1.4826*MAD
        reject.epoch.zValEpochHFNoise(ichan,:) = (noisiness - noisinessMedian) ./ noisinessSD;
    end
    reject.epoch.epochHFNoise = reject.epoch.zValEpochHFNoise > 2 * reject.HFAnalysis.threshold;
    
    
    
    
    %% METHOD 2.3: High Frequency (HF) Noise detection: across channels; all time points (X is transformed to 2D form: Nchannel * (Ntp*Nepoch)
    noisiness = mad(reshape(X,Nchan,[])- reshape(Xlp,Nchan,[]), 1, 2)./mad(reshape(Xlp,Nchan,[]), 1, 2);
    noisinessMedian = nanmedian(noisiness);
    noisinessSD = mad(noisiness, 1, 1)*1.4826; % for a normal distribution sigma=1.4826*MAD
    reject.chan.zValHFNoise = (noisiness - noisinessMedian) ./ noisinessSD;
    reject.chan.HFNoise = (reject.chan.zValHFNoise > reject.HFAnalysis.threshold);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %% Final Bad Channel detection
    
    chanReject1 = reject.chan.HFNoise;
    chanReject2 = reject.chan.chanDeviation;
    chanReject3 = (sum(reject.epoch.chanHFNoise, 2)/Nepoch > 0.2);
    chanReject4 = sum(reject.epoch.chanDeviation,2)/Nepoch > 0.2;
    reject.chanReject = [chanReject1 chanReject2 chanReject3 chanReject4];
    
    if par.HFAnalysis.channelRejection
    reject.chanRejectFinal = chanReject1| chanReject2 | chanReject3 | chanReject4;
    else
    reject.chanRejectFinal = chanReject2| chanReject4;    
    end
    
    %%
    
end


%%% Removing detected bad channels for epoch detection
XX = X(~reject.chanRejectFinal, :, :);
Nchan2 = size(XX,1);






%% %% BAD EPOCH DETECTION


if any(reject.chanRejectFinal)
    
    if par.deviationAnalysis.run
        %% METHOD 1.1:  Deviation across channels; epoch by epoch
        reject.epoch.zValChanDeviation2 = zeros(Nchan2, Nepoch);
        for iepoch = 1:Nepoch
            chanDeviation = 0.7413 * iqr(squeeze(XX(:,:,iepoch)), 2); % Robust estimate of SD
            chanDeviationSD =  0.7413 * iqr(chanDeviation);
            chanDeviationMedian = nanmedian(chanDeviation);
            reject.epoch.zValChanDeviation2(:,iepoch) = ...
                (chanDeviation - chanDeviationMedian) / chanDeviationSD;
        end
        reject.epoch.chanDeviation2 = (abs(reject.epoch.zValChanDeviation2) > reject.deviationAnalysis.threshold);
        reject.epoch.chanExDeviation2 = (abs(reject.epoch.zValChanDeviation2) > 2*reject.deviationAnalysis.threshold);
        
        
        %% METHOD 1.2:  Deviation across epochs; channel by channel
        reject.epoch.zValEpochDeviation2 = zeros(Nchan2, Nepoch);
        for ichan = 1:Nchan2
            epochDeviation = 0.7413 * iqr(squeeze(XX(ichan,:,:)), 1); % Robust estimate of SD
            epochDeviationSD =  0.7413 * iqr(epochDeviation);
            epochDeviationMedian = nanmedian(epochDeviation);
            reject.epoch.zValEpochDeviation2(ichan,:) = ...
                (epochDeviation - epochDeviationMedian) / epochDeviationSD;
        end
        reject.epoch.epochDeviation2 = (abs(reject.epoch.zValEpochDeviation2) > reject.deviationAnalysis.threshold);
        reject.epoch.epochExDeviation2 = (abs(reject.epoch.zValEpochDeviation2) > 2 * reject.deviationAnalysis.threshold);
        
    end
    
    
    if par.HFAnalysis.run
        
        XXlp = Xlp(~reject.chanRejectFinal, :, :);
        
        %% METHOD 2.1:  High Frequency (HF) Noise detection: across channels; epoch by epoch
        reject.epoch.zValHFNoise2 = zeros(Nchan2, Nepoch);
        for iepoch = 1:Nepoch
            noisiness = mad(XX(:,:,iepoch)- XXlp(:,:,iepoch), 1, 2)./mad(XXlp(:,:,iepoch), 1, 2);
            noisinessMedian = nanmedian(noisiness);
            noisinessSD = mad(noisiness, 1, 1)*1.4826; % for a normal distribution sigma=1.4826*MAD
            reject.epoch.zValChanHFNoise2(:,iepoch) = (noisiness - noisinessMedian) ./ noisinessSD;
        end
        reject.epoch.chanHFNoise2 = reject.epoch.zValChanHFNoise2 > reject.HFAnalysis.threshold;
        
        
        %% METHOD 2.2:  High Frequency (HF) Noise detection: across channels; epoch by epoch
        reject.epoch.zValEpochHFNoise2 = zeros(Nchan2, Nepoch);
        for ichan = 1:Nchan2
            noisiness = mad(squeeze(XX(ichan,:,:)- XXlp(ichan,:,:)), 1, 1)./mad(squeeze(XXlp(ichan,:,:)), 1, 1);
            noisinessMedian = nanmedian(noisiness);
            noisinessSD = mad(noisiness, 1, 2)*1.4826; % for a normal distribution sigma=1.4826*MAD
            reject.epoch.zValEpochHFNoise2(ichan,:) = (noisiness - noisinessMedian) ./ noisinessSD;
        end
        reject.epoch.epochHFNoise2 = reject.epoch.zValEpochHFNoise2 > reject.HFAnalysis.threshold;
    end
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %% Final: Bad Epochs detection
if any(reject.chanRejectFinal)
    epochReject1 = logical(sum((and(reject.epoch.chanHFNoise2,reject.epoch.chanDeviation2))));
    epochReject2 = sum(reject.epoch.epochHFNoise2, 1)/Nchan2 > 0.1;
    epochReject3 = sum(reject.epoch.epochDeviation2, 1)/Nchan2 > 0.1;
    epochReject4 = logical(sum(reject.epoch.epochExDeviation2));
    epochReject5 = logical(sum(reject.epoch.chanExDeviation2));
    reject.epochReject = [epochReject1;epochReject2;epochReject3;epochReject4;epochReject5];
    reject.epochRejectFinal = epochReject1 | epochReject2 | epochReject3 | (epochReject4 & epochReject5);
else
    epochReject1 = logical(sum((and(reject.epoch.chanHFNoise,reject.epoch.chanDeviation))));
    epochReject2 = sum(reject.epoch.epochHFNoise, 1)/Nchan > 0.1;
    epochReject3 = sum(reject.epoch.epochDeviation, 1)/Nchan > 0.1;
    epochReject4 = logical(sum(reject.epoch.epochExDeviation));
    epochReject5 = logical(sum(reject.epoch.chanExDeviation));
    reject.epochReject = [epochReject1;epochReject2;epochReject3;epochReject4;epochReject5];
    reject.epochRejectFinal = epochReject1 | epochReject2 | epochReject3 | (epochReject4 & epochReject5);
end



