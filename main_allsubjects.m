% Calls main.m for all subjects. Refer to that function for documentation.

function main_allsubjects(nshuf, varargin)
    
    % add path
    addpath_cluster

    eeglab
    g = finputcheck(varargin, { ...
        'n'              'integer'       { }                5;
        'alpha'          'float'         { }              0.05;
        'freq_manual'    'string'        { 'on' 'off' }   'off';
        'f1'             'integer'       { }              11; 
        'f2'             'integer'       { }              22; 
        'run_ica'        'string'        { 'on' 'off' }   'off';
        'poolsize'       'integer'       { }              1;
        });
    if ischar(g), error(g); end

    subjects = 1:1:37;
    exclude = [1, 2, 6, 7, 10, 13, 20, 24, 26, 32, 36]; % 26 subjects
    
    for isub = 1:length(subjects)
        tic
        if ismember(isub, exclude)
            continue
        else
            fprintf('Subject %d  .............................................................................................. \n', isub)
            main(nshuf, isub, 'n', g.n, 'alpha', g.alpha, 'freq_manual', g.freq_manual, 'f1', g.f1, 'f2', g.f2, 'run_ica', g.run_ica, 'poolsize', g.poolsize)
        end
        toc
    end
end