% Script to preprocess VitalBCI motor imagery data. Script is based on https://github.com/fpellegrini/MotorImag/blob/main/fp_preprocess.m 
% by Franziska Pellegrini.

%% Add paths
addpath('/Users/nguyentiendung/Desktop/Studium/Charite/matlab/eeglab/')
eeglab
addpath('/Users/nguyentiendung/Desktop/Studium/Charite/matlab/eeglab/plugins/dipfit/standard_BEM/elec/')

%% Set up directories and load classification scores
DIRIN = '/Volumes/Seagate/VitalBCI/mat/imag/';
DIROUT = '/Volumes/Seagate/MotorImag/data/';
load('/Volumes/Seagate/VitalBCI/mat/scores.mat')

%% Run preprocessing pipeline
nsub = 40;
DB_noise = db_motorImag_noisechans;

for isub = 2:nsub % subject 1 is not available
    sub = ['vp' num2str(isub)];
    if strcmp(scores{isub}(1),'1') % if subject has good classification performance
        % load data
        load([DIRIN sub '.mat'])
        data = 0.1 * double(cnt.x)';
    
        % filter before resampling
        fs = 1000;
        notchb = [48 52];
        [bband, aband] = butter(2, notchb/fs*2, 'stop');
        [bhigh, ahigh] = butter(2, 1/fs*2, 'high');
        [blow, alow] = butter(2, 45/fs*2, 'low');
        
        data = filtfilt(bhigh, ahigh, data');
        data = filtfilt(bband, aband, data);
        data = filtfilt(blow, alow, data');
        
        % resample to 100 Hz
        dsRatio = 10;
        data = data(:, 1:dsRatio:end);
        fs = 100;
        
        % enter marker positions
        if ~(strcmp(mrk.className{1},'left')) || ~(strcmp(mrk.className{2},'right'))
            error(['check classNames in subject ' subs{isub}])
        end
        data(end+1, :) = zeros(1,size(data,2));
        data(end, mrk.pos) = mrk.toe;
        
        % load data into EEG struct
        EEG = pop_importdata('dataformat', 'array', 'nbchan', 0, 'data', 'data', 'setname', 'test', 'srate', fs, 'pnts', 0, 'xmin', 0);
        EEG = eeg_checkset(EEG);
        
        % add events
        EEG = pop_chanevent(EEG, EEG.nbchan, 'edge', 'leading');
        
        % load template to add channel information
        eloc = readlocs('standard_1005.elc');
        
        % remove all channels of eloc that are not contained in the data
        out = [];
        out_eloc = [];
        in = [];
        for ii = 1:length(eloc)
            if ~ismember(eloc(ii).labels,cnt.clab)
                out = [out,ii];
            else
                in = [in find(strcmp(eloc(ii).labels,cnt.clab))];
            end
        end
        eloc(out)=[];
        
        % remove channels in data that have not been found in eloc template
        out_labs = setdiff(1:EEG.nbchan, in);
        EEG = pop_select(EEG, 'nochannel', out_labs);
        
        % bring channels in right order
        eloc1 = eloc;
        for ii = 1:length(eloc)
            name_eloc = eloc(ii).labels;
            ind_clab = find(strcmp(eloc(ii).labels,cnt.clab));
            eloc1(ind_clab) = eloc(ii);
        end
        
        % add chanlocs to EEG struct
        EEG = pop_editset(EEG, 'chanlocs', eloc1, 'setname', 'setName');
        EEG = eeg_checkset(EEG);
        
        eeglabp = fileparts(which('eeglab.m'));
        EEG = pop_chanedit(EEG, 'lookup', fullfile(eeglabp, 'plugins','dipfit','standard_BEM','elec','standard_1005.elc'));  
        chanlocs = EEG(1).chanlocs;
        
        % manual rejection of bad channels 
        clear chans
        chans = matchdbs_motorImag(DB_noise, sub);
        if ~isempty(chans)
            EEG = pop_select(EEG, 'nochannel', chans);
            EEG = pop_interp(EEG, chanlocs);
        end
        
        % save preprocessed EEG struct
        pop_saveset(EEG,'filename',['prep_' sub],'filepath',DIROUT)
    end
end

