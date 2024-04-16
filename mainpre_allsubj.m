% Calls main_preanalysis.m for selected subjects with their individual f1. Refer to that function for documentation.

function mainpre_allsubj(n_shuf, varargin)
    % add path
    addpath_cluster

    eeglab
    g = finputcheck(varargin, { ...
        'alpha'          'float'         { }              0.05;
        'poolsize'       'integer'       { }              1;
        'epleng'         'integer'       { }              2;
        'downsample'     'string'        {'on' 'off'}     'on';
        'freq_down'      'integer'       { }              125;
        });
    if ischar(g), error(g); end
    
    subjects = [317, 511, 429, 306, 355, 360, 441, 330, 385, 413];
    f1 = [9.6, 9.4, 9.4, 9.1, 9.8, 10.4, 10.4, 9.9, 9.4, 9.2];

    for ind_sub = 1:length(subjects)
        g.f1 = f1(ind_sub);
        main_preanalysis(n_shuf, subjects(ind_sub), g)
    end
end
