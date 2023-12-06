
function [hdr,dat] = wlb_read_nextim(filename)

[hdr] = read_nextim_nxh(filename);

%[dat] = read_nextim_nxe(filename, hdr, begsample, endsample, chanindx);
[dat] = read_nextim_nxe(filename, hdr);
end

function [dat] = read_nextim_nxe(filename, hdr)

% FIXME it would be nice to also implement the efficient reading of the
% selected channels for the other file formats but at the moment only the
% implementation of the binary multiplexed formats is smart enough.
begsample = 1;
endsample = hdr.nSamples;

switch lower(hdr.BinaryFormat)
    case 'int_16'
        sampletype = 'int16';
        samplesize = 2;
    case 'int_32'
        sampletype = 'int32';
        samplesize = 4;
    case 'ieee_float_32'
        sampletype = 'float32';
        samplesize = 4;
end % case

fid = fopen(filename, 'rb', 'ieee-le');

% read all the channels
fseek(fid, hdr.NumberOfChannels*samplesize*(begsample-1), 'cof');
dat = fread(fid, [hdr.NumberOfChannels, (endsample-begsample+1)], sampletype);
% compute real microvolts using the calibration factor (resolution)
% calib = diag(hdr.resolution);
% % using a sparse multiplication speeds it up
% dat = full(sparse(calib) * dat);
calib = reshape(hdr.resolution,[],1);
for k = 1:size(dat,2)
    dat(:,k) = calib.*dat(:,k);
end  

fclose(fid);

end

function [hdr] = read_nextim_nxh(filename)

[p, f, x] = fileparts(filename);

hdr.DataFile         = [f,x];
hdr.MarkerFile       = [];
hdr.DataFormat       = 'BINARY';
hdr.DataOrientation  = 'MULTIPLEXED';
hdr.BinaryFormat     = 'INT_16';
hdr.NumberOfChannels = 64;
hdr.Fs = 1450; %Hz
hdr.SamplingInterval = 1e6/hdr.Fs;   % microseconds

hdr.label = {'TRG';'TRG';'TRG';'EOG';'FP1';'FPz';'FP2';'AF1';'AFz';'AF2';'F7';...
    'F5';'F1';'Fz';'F2';'F6';'F8';'FT9';'FT7';'FC5';'FC3';'FC1';'FCz';'FC2';'FC4';...
    'FC6';'FT8';'FT10';'T3';'C5';'C3';'C1';'Cz';'C2';'C4';'C6';'T4';'TP9';'TP7';...
    'CP5';'CP3';'CP1';'CPz';'Cp2';'Cp4';'CP6';'TP8';'TP10';'P9';'P7';'P3';'P1';...
    'Pz';'P2';'P4';'P8';'P10';'PO3';'POz';'PO4';'O1';'Oz';'O2';'Iz'};

hdr.reference = cell(size(hdr.label,1),size(hdr.label,2));

sfEEG = (1/2000) * (10/65535) * 1000000;
sfEOG = (1/400) * (10/65535) * 1000000;
sfTRIG = 10 * (10/65535);

hdr.resolution = [sfTRIG ;sfTRIG ;sfTRIG ;sfEOG ;repmat(sfEEG,60,1)];

hdr.chanunit = repmat({'µV'}, size(hdr.label));

hdr.chantype = ['trg';'trg';'trg';'eog';repmat({'eeg'}, 60,1)];
    
sph_theta_besa = [-130;-130;-130;-130;-92;92;92;-74;69;74;-92;-75;-50;46;50;75;92;...
    -115;-92;-72;-50;-32;23;32;50;72;92;115;-92;-69;-46;-23;0;23;46;69;92;-115;-92;...
    -72;-50;-32;23;32;50;72;92;115;-115;-92;-60;-50;46;50;60;92;115;-74;69;74;-92;92;...
    92;115];
sph_phi_besa = [-45;-45;-45;-45;-72;90;72;-65;90;65;-36;-41;-68;90;68;41;36;-15;-18;...
    -21;-28;-45;90;45;28;21;18;15;0;0;0;0;0;0;0;0;0;15;18;21;28;45;-90;-45;-28;-21;...
    -18;-15;36;36;51;68;-90;-68;-51;-36;-36;65;-90;-65;72;-90;-72;-90];

for i=1:hdr.NumberOfChannels
    hdr.layout.pos(i).type = hdr.chantype(i);
    hdr.layout.pos(i).labels = hdr.label(i);
    hdr.layout.pos(i).sph_theta_besa = sph_theta_besa(i);
    hdr.layout.pos(i).sph_phi_besa = sph_phi_besa(i);
end

% determine the number of samples by looking at the binary file
Datafile = fullfile(p, hdr.DataFile); % add full-path to datafile

% the data file is supposed to be located in the same directory as the header file
% but that might be on another location than the present working directory
info = dir(Datafile);
if isempty(info)
    error('cannot determine the location of the data file %s', hdr.DataFile);
end
switch lower(hdr.BinaryFormat)
    case 'int_16';
        hdr.nSamples = info.bytes ./ (hdr.NumberOfChannels*2);
    case 'int_32';
        hdr.nSamples = info.bytes ./ (hdr.NumberOfChannels*4);
    case 'ieee_float_32';
        hdr.nSamples = info.bytes ./ (hdr.NumberOfChannels*4);
end

if isinf(hdr.nSamples)
    warning('cannot determine number of samples for this sub-fileformat');
end

% the number of trials is unkown, assume continuous data
hdr.nTrials     = 1;
hdr.nSamplesPre = 0;

% ensure that the labels are in a column
hdr.label      = hdr.label(:);
hdr.reference  = hdr.reference(:);
hdr.resolution = hdr.resolution(:);

end
