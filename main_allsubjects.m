% Calls main.m for all subjects. Refer to that function for documentation.

function main_allsubjects(nshuf, varargin)
    
    % add path
    addpath_cluster

    eeglab
    g = finputcheck(varargin, { ...
        'n'              'integer'       { }                5;
        'alpha'          'float'         { }              0.05;
        'freq_manual'    'string'        { 'on' 'off' }   'on';
        'run_ica'        'string'        { 'on' 'off' }   'off';
        'poolsize'       'integer'       { }              1;
        });
    if ischar(g), error(g); end

    subjects = [317, 511, 429, 306, 355, 360, 441, 330, 385, 413];
    f1 = [9.6, 9.4, 9.4, 9.1, 9.8, 10.4, 10.4, 9.9, 9.4, 9.2];

    for isub = 1:length(subjects)
        disp(['Subject ' int2str(isub) '..............................................................................................'])
        % (f1, f1, 2*f1)
        tic
        main(nshuf, subjects(isub), 'n', g.n, 'alpha', g.alpha, 'freq_manual', g.freq_manual, 'f1', f1(isub), 'f2', f1(isub), 'run_ica', g.run_ica, 'poolsize', g.poolsize, 'antisymm', [3 2 1])
        toc

        % (f1, f2, f1+f2)
        tic
        main(nshuf, subjects(isub), 'n', g.n, 'alpha', g.alpha, 'freq_manual', g.freq_manual, 'f1', f1(isub), 'f2', 2 * f1(isub), 'run_ica', g.run_ica, 'poolsize', g.poolsize, 'antisymm', [2 1 3])
        toc
    end
end