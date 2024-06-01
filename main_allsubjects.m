% Calls main.m for all subjects. Refer to that function for documentation.

function main_allsubjects(nshuf, varargin)
    
    % add path
    addpath_cluster

    eeglab
    g = finputcheck(varargin, { ...
        'n'                'integer'       { }                5;
        'alpha'            'float'         { }              0.05;
        'freq_manual'      'string'        { 'on' 'off' }   'on';
        'run_ica'          'string'        { 'on' 'off' }   'off';
        'poolsize'         'integer'       { }              1;
        'epleng'           'integer'       { }              2;
        'downsample'       'string'        { 'on' 'off'}    'on';
        'freq_down'        'integer'       { }              125;
        'bispec_type'      'string'        { }              ''; 
        'antisymm'         'integer'       { }              [1 2 3];
        'train_test'       'string'        { 'on' 'off'}    'off';
        'dim_chan'         'integer'       { 1 2 3 }         2;
        });
    if ischar(g), error(g); end

    subjects = [317, 511, 429, 306, 355, 360, 441, 330, 385, 413];
    f1 = [9.6, 9.4, 9.4, 9.1, 9.8, 10.4, 10.4, 9.9, 9.4, 9.2];
    % subjects = [511, 429, 306, 355, 360, 441, 330, 385, 413];
    % f1 = [9.4, 9.4, 9.1, 9.8, 10.4, 10.4, 9.9, 9.4, 9.2];

    for isub = 1:length(subjects)
        disp(['Subject ' int2str(isub) '..............................................................................................'])
        g.f1 = f1(isub);
        g.f2 = g.f1;
        tic
        main(nshuf, subjects(isub), g)
        toc
    end
end