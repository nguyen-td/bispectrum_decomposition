% Low-level function to compute power spectra using EEGLAB's function. The title will have the form 
% "condition subject", so you can pick any combination that works.
% 
% Inputs:
%   EEG        - EEG struct
%   condition  - condition for title, e.g., eyes open (EO) or eyes closed (EC), can be any other label too 
%   subject    - subject ID for title, can be any other label too
%   result_dir - result directory
%
% Optional inputs:
%   title_str  - manual title, if not passed the default combination "subject condition" is being used. Do not add the .png extension here.
%   f1         - [integer] Frequency for plotting topomaps. The function will then plot its first and second harmonic respectively.

function plot_spectra(EEG, condition, subject, result_dir, varargin)

    g = finputcheck(varargin, { ...
        'title_str'      'string'        []              [];
        'f1'             'integer'       { }             10;
        });
    if ischar(g), error(g); end

    EEG = eeg_checkset( EEG );
    f = figure; 
    pop_spectopo(EEG, 1, [0 EEG.xmax * 1000], 'EEG' , 'percent', 80, 'freq', [g.f1 g.f1*2 g.f1*3], 'freqrange',[0 EEG.srate/2], ...
        'electrodes', 'on', 'overlap', 50, 'title', [condition, ' ', subject]); 

    if isempty(g.title_str)
        exportgraphics(gcf, [result_dir, subject, '_', condition, '_psd.png']);
    else
        exportgraphics(gcf, [result_dir, g.title_str '.png']);
    end
    close(f)

end