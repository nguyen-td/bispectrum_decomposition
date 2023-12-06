%% Add paths
addpath('/Users/nguyentiendung/Desktop/Studium/Charite/matlab/eeglab') % eeglab path
addpath(genpath('/Volumes/PortableSSD/EmergencyBreaking')) % data location
addpath(genpath('/Users/nguyentiendung/Desktop/Studium/Charite/Research/Project 1/bispectrum_decomposition/EmergencyBreaking')) % script location
addpath('/Users/nguyentiendung/Desktop/Studium/Charite/Matlab/eeglab/plugins/Fieldtrip-lite20230125/template/electrode')

% data directory
data_dir = '/Volumes/PortableSSD/EmergencyBreaking/data/';

%%% Creating result folder (Plots will be saved here)
resultDir = ['analysis_output/preprocessing/'];
if ~isdir(resultDir)
    mkdir(resultDir)
end


%%% Creating output data (matfile) folder
outputDataDir = [resultDir, 'data/'];
if ~isdir(outputDataDir)
    mkdir(outputDataDir)
end
    
load metadata

fs = 100;
fres = fs;

frqs = sfreqs(fres, fs);
% maxfreq_ind = max(find(frqs <= 45));
% frqs = frqs(1:maxfreq_ind);

%% Preprocessing Settings
saveOutputData = 1;   % 1: To save preprocessd data and preprocessing details  0: Not save

%% %%%%%%%%%%%%%%%%%
lpFilter =   45;       % low-pass filter cut-off
hpFilter =   1;      % high-pass filter cut-off
bsFilter =   [49 51];       % band-stop filter range
filterOrder = 2;    % Butterworth filter order
dsRatio =  2;       % downsampling rate
%dsRatio = fs/EGE.rate
epochLeng = 2;   % epoch length in seconds
mincomp = 30;

% Updating 'info' variable
c = fix(clock);
info = [];
info.prepStartTime = [' ',num2str(c(3)),'.',num2str(c(2)),'.',num2str(c(1)),'     ',num2str(c(4)),':',num2str(c(5))];
info.prep.lpFilter = lpFilter;
info.prep.hpFilter  = hpFilter;
info.prep.filterOrder = filterOrder;


[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

nsub = 18;
exclude = [15]; 

isub = 1; % for testing
for isub = 1:nsub
    if ismember(isub, exclude)
        continue
    end

    fprintf('Current subject: %d\n', isub)
    sub = ['vp' num2str(isub)];
    load([data_dir sub '.mat'])
    data = 0.1 * double(cnt.x)'; % muV?
    
    % remove non-EEG channels
    chans_reject = ["EOGv", "EOGh", "EMGf", "lead_gas", "lead_brake", "dist_to_lead", "wheel_X", "wheel_Y", "gas", "brake"];
    idx_reject = find(ismember(cnt.clab, chans_reject));
    data(idx_reject, :) = [];
    cnt.clab(idx_reject) = [];

    % enter marker positions
    data(end+1, :) = zeros(1, size(data, 2));
    mrk_pos = int64((mrk.time / 1000) * fs); % convert mrk.time from ms to seconds and then to sampling points
    data(end, mrk_pos) = onehotdecode(mrk.y, [1, 2, 3, 4, 5], 1);

    % load data into EEG struct
    EEG = pop_importdata('dataformat', 'array', 'nbchan', 0, 'data', 'data', 'setname', 'test', 'srate', cnt.fs, 'pnts', 0, 'xmin', 0);
    EEG = eeg_checkset(EEG);
    EEG.epoch = 1;

    % add events
    EEG = pop_chanevent(EEG, EEG.nbchan, 'edge', 'leading');

    % load template to add channel information
    eloc = readlocs('standard_1005.elc');
    
    % remove all channels of eloc that are not contained in the data
    out = [];
    out_eloc = [];
    in = [];
    for ii = 1:length(eloc)
        if ~ismember(eloc(ii).labels, cnt.clab)
            out = [out, ii];
        else
            in = [in find(strcmp(eloc(ii).labels, cnt.clab))];
        end
    end
    eloc(out)=[];
    
    % bring channels in right order
    eloc1 = eloc;
    for ii = 1:length(eloc)
        name_eloc = eloc(ii).labels;
        ind_clab = find(strcmp(eloc(ii).labels, cnt.clab));
        eloc1(ind_clab) = eloc(ii);
    end
    
    % add chanlocs to EEG struct
    EEG = pop_editset(EEG, 'chanlocs', eloc1, 'setname', 'setName');
    EEG = eeg_checkset(EEG);

    % Updating "info" variable
    info.prep.Nchan0 = EEG.nbchan; info.data.fs0 = EEG.srate;
    info.data.lengthPoints = EEG.pnts;
    info.prep.ref0 = 'common'; info.prep.chanLabels = {EEG.chanlocs.labels};

    %% Plotting Channel locations
    EEG = eeg_checkset( EEG );
    figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
    saveas(gcf,[resultDir,'topo_labels'],'jpg');
    EEG = eeg_checkset( EEG );
    figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'numpoint', 'chaninfo', EEG.chaninfo);
    saveas(gcf,[resultDir,'topo_indices'],'jpg');

    %% Plotting Spectra of raw data
    EEG = eeg_checkset( EEG );
    figure; pop_spectopo(EEG, 1, [0 EEG.xmax*1000], 'EEG' , 'percent', 15, 'freq', [10 20], 'freqrange',[0 45],'electrodes','on');
    saveas(gcf,[resultDir,'spectra_raw'],'jpg');

    %% zero-phase filtering
    [b_low, a_low] = butter(filterOrder, lpFilter/(EEG.srate/2), 'low');
    [b_high, a_high] = butter(filterOrder, hpFilter/(EEG.srate/2), 'high');
    [b_notch, a_notch] = butter(filterOrder, bsFilter/(EEG.srate/2),'stop');
    a_all = poly([roots(a_low);roots(a_high); roots(a_notch)]);
    b_all = conv(conv(b_low,b_high), b_notch);
    %     fvtool(b_all, a_all)
    EEG.data = filtfilt(b_all, a_all, double(EEG.data)')';

    %% plain downsampling to 100 Hz. 
    EEG.data = EEG.data(:, 1:dsRatio:end);
    EEG.srate = EEG.srate/dsRatio;
    EEG.pnts    = size(EEG.data,2);
    EEG.xmax    = EEG.xmin + (EEG.pnts-1)/EEG.srate;
    EEG.times   = linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts);
    EEG = eeg_checkset( EEG );

    %% Epoching: Trasforming data from Nchan*NtpAll  to Nchan*NtpEpoch*Nepoch
    %epochLeng = 2; % epoch length in seconds
    [Nchan, Ntp]=size(EEG.data);
    EEG = pop_select( EEG,'point',[1 (Ntp - mod(Ntp,EEG.srate*epochLeng))] );
    Nepoch = EEG.pnts / (EEG.srate*epochLeng);
    EEG = eeg_checkset(EEG);
    %EEG.data = EEG.data(:, 1:(Ntp - mod(Ntp,fsNew*8)) );
    for ievent = 1: Nepoch
        EEG.event(ievent).type = num2str(epochLeng);
        EEG.event(ievent).latency = (ievent-1)*(EEG.srate*epochLeng)+1;
        EEG.event(ievent).duration = epochLeng;
    end
    EEG = eeg_checkset( EEG );
    EEG = pop_epoch( EEG, {  num2str(epochLeng)  }, [0  epochLeng], 'epochinfo', 'yes');

    %% Plotting Spectra of filtered and downsampled data
    EEG = eeg_checkset( EEG );
    figure; pop_spectopo(EEG, 1, [0      EEG.xmax*1000*Nepoch], 'EEG' , 'percent', 15, 'freq', [10 20 ], 'freqrange',[0 45],'electrodes','on');
    saveas(gcf,[resultDir,'spectra_prep1'],'jpg');

    % Updating 'info' variable
    % info.subject.ID = isbj;
    % info.subject.eyeCondition = EYE;
    % info.subject.medicalCondition = MED;
    % info.subject.fileDirectory = dataDir;
    info.subject.outputDirectory = resultDir;
    info.prep.fs = EEG.srate;  % sampling frequency after down-sampling
    info.prep.Nepochs0 = size(EEG.data,3);
    info.prep.epochPoints = size(EEG.data, 2);
    info.prep.Nchannels0 = size(EEG.data, 1);
    info.prep.dsRation = dsRatio;
    info.prep.ref = '';

    %% Detecting Bad channels/epochs
    par.HFAnalysis.run = 1 ;
    par.HFAnalysis.channelRejection = 1;
    info.outlier = k1_detect_bad_epoch_channel(EEG,par);

    %% Detecting channels/epochs rejected because of Strong alpha activity
    [b, a] = butter(2, [5 12]/(EEG.srate/2),'stop');
    XNoAlpha = reshape(filtfilt(b, a, reshape(double(EEG.data), Nchan, [])')',Nchan,[],Nepoch);
    par = [];
    par.srate = EEG.srate;
    par.deviationAnalysis.run = 1;
    par.HFAnalysis.run = 1;
    par.deviationAnalysis.threshold = 4.5;
    par.HFAnalysis.channelRejection =1;
    infoAlpha = k1_detect_bad_epoch_channel(XNoAlpha, par);

    info.outlier.epochRejectFinal = info.outlier.epochRejectFinal & infoAlpha.epochRejectFinal;  

    %% Plotting Rejceted channels/epochs
    epoch = find(info.outlier.epochRejectFinal);
    epochchanind = cell(1, length(epoch)); %{[],[3,9,12],[6:12,16,18]};
    rejepochcol =  [.6 .8 .9];
    rejepoch = zeros(1, EEG.trials);
    rejepoch(epoch) = ones(1, length(epoch));
    rejepochE = zeros(EEG.nbchan, EEG.trials);
    for i = 1:length(find(rejepoch))
      rejepochE(epochchanind{i},epoch(i))=ones(size(epochchanind{i}));
    end
    winrej = trial2eegplot(rejepoch, rejepochE, EEG.pnts, rejepochcol);
    colors = cell(1,Nchan); colors(1,:) = {'k'};
    colors(1,info.outlier.chanRejectFinal) = {'r'};
    eegplot(EEG.data, 'srate', EEG.srate, 'title', 'Rejected Channels Epochs', ...
      'limits', [EEG.xmin EEG.xmax]*1000, 'color', colors, 'winrej',winrej, 'eloc_file', EEG.chanlocs);

    %% Additional plots for testing channel/epoch rejection
    if 1
      figure, imagesc(info.outlier.epoch.epochDeviation);
      title('Epoch-Epoch-Deviation')
      xlabel('Epoch','fontweight','b'); ylabel('Channel','fontweight','b');
    %         xticks(5:5:Nepoch);
    %         yticks(1:2:Nchan);
      saveas(gcf,[resultDir,'epoch_epoch_deviation'],'jpg');
    
      figure, imagesc(info.outlier.epoch.chanDeviation);
      title('Epoch-Chan-Deviation')
      xlabel('Epoch','fontweight','b'); ylabel('Channel','fontweight','b');
    %         xticks(5:5:Nepoch);
    %         yticks(1:2:Nchan);
      saveas(gcf,[resultDir,'epoch_chan_deviation'],'jpg');
    
      figure, imagesc(info.outlier.epoch.epochHFNoise);
      title('Epoch-Epoch-HFNOISE')
      xlabel('Epoch','fontweight','b'); ylabel('Channel','fontweight','b');
    %         xticks(5:5:Nepoch);
    %         yticks(1:2:Nchan);
      saveas(gcf,[resultDir,'epoch_epoch_hfnoise'],'jpg');
    
      figure, imagesc(info.outlier.epoch.chanHFNoise);
      title('Epoch-Chan-HFNOISE')
      xlabel('Epoch','fontweight','b'); ylabel('Channel','fontweight','b');
    %         xticks(5:5:Nepoch);
    %         yticks(1:2:Nchan);
      saveas(gcf,[resultDir,'epoch_chan_hfnoise'],'jpg');
    
      chanEpochMatrix = zeros(Nchan, Nepoch);
      chanEpochMatrix(info.outlier.chanRejectFinal, :) = 1;
      chanEpochMatrix(:, info.outlier.epochRejectFinal) = 1;
      figure, imagesc(chanEpochMatrix);
      title('Detected outlier channels epochs')
      xlabel('Epoch','fontweight','b'); ylabel('Channel','fontweight','b');
    %         xticks(5:5:Nepoch);
    %         yticks(1:2:Nchan);
      saveas(gcf,[resultDir,'outlier_chans_epochs'],'jpg');
    end

    %% Rejecting bad channels/epochs before ICA
    info.originalChanlocs = EEG.chanlocs;
    EEG = pop_select( EEG,'notrial',find(info.outlier.epochRejectFinal) ,'nochannel',find(info.outlier.chanRejectFinal));

    %% Plotting Spectra after channel/epochs rejection
    EEG = eeg_checkset( EEG );
    figure; pop_spectopo(EEG, 1, [0      EEG.xmax*1000*Nepoch], 'EEG' , 'percent', 15, 'freq', [10 20 ], 'freqrange',[0 45],'electrodes','on');
    saveas(gcf,[resultDir,'spectra_prep2'],'jpg');
    
end