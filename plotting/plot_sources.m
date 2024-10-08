% Plot and save mixed and demixed sources as brain plots.
%
% Inputs:
%   F_moca      - (n_voxels x n_dum x n) demixed sources
%   n           - number of estimated sources 
%   cortex75k   - 75k cortex structure from the NYhead for later plotting. Will be ignored if you use a Brainstorm cortex.
%   cortex2k    - 2k subset of the cortex structure from the NYhead for later plotting. Will be ignored if you use a Brainstorm cortex.
%   roi_inds    - (n x 1) array of ROI indicess
%   cm          - colormap
%   isub        - subject ID, used for cortex plot
%   DIROUT      - output directory to save images
%
% Optional input:
%   bispec_type       - [string] type of bispectrum (for file name), default is '' (empty string)
%   cortex_BS         - [struct] Brainstorm cortex. If passed, the cortex75k and cortex2k will be ignored.
%   in_normal_to_high - (1 x 15k) array inverse vector to plot sources on the high-resolution cortex. Only for Brainstorm-based cortices.

function plot_sources(F_moca, n, cortex75k, cortex2k, roi_inds, cm, isub, DIROUT, varargin)

    g = finputcheck(varargin, { ...
        'bispec_type'        'string'     { }             '';
        'cortex_BS'          'struct'     { }             struct;
        'in_normal_to_high'  'float'      { }             [];
        'allseeds'           'string'     { 'on' 'off' }  'off';
         });
    if ischar(g), error(g); end

    % get max. absolute value
    max_val = max(sum(F_moca.^2, 2), [], 'all');
    
    % plot and save demixed source
    close all
    for i = 1:n
        source_moca = sum(F_moca(:, :, i).^2, 2); 
        f_name_moca = [DIROUT 'F' g.bispec_type '_demixed' int2str(i) '_' int2str(isub) '_'];
        if isempty(fieldnames(g.cortex_BS)) % use New York head
            if strcmpi(g.allseeds, 'off')
                try
                    allplots_cortex_nyhead_v2(cortex75k, source_moca(cortex2k.in_to_cortex75K_geod), [0 max_val], cm, 'a.u.', 0, f_name_moca, ...
                        {roi_inds(i, :)});
                catch
                    allplots_cortex_nyhead_v2(cortex75k, source_moca(cortex2k.in_to_cortex75K_geod), [0 max_val], cm, 'a.u.', 0, f_name_moca);
                end
            else
                allplots_cortex_nyhead_v2(cortex75k, source_moca(cortex2k.in_to_cortex75K_geod), [0 max_val], cm, 'a.u.', 0, f_name_moca, ...
                        {roi_inds});
            end
        else % use Brainstorm-based head
            if isempty(g.in_normal_to_high)
                error('The "in_normal_to_high" is missing for plotting sources on the high-resolution cortex.')
            else
                if strcmpi(g.allseeds, 'off')
                    try
                        allplots_cortex_BS_v2(g.cortex_BS, source_moca(g.in_normal_to_high), [0 max_val], cm, 'a.u.', 0.1, f_name_moca, ...
                            {roi_inds(i, :)});
                    catch
                        allplots_cortex_BS_v2(g.cortex_BS, source_moca(g.in_normal_to_high), [0 max_val], cm, 'a.u.', 0.1, f_name_moca);
                    end
                else
                    allplots_cortex_BS_v2(g.cortex_BS, source_moca(g.in_normal_to_high), [0 max_val], cm, 'a.u.', 0.1, f_name_moca, ...
                            {roi_inds});
                end
            end
        end
    end

end