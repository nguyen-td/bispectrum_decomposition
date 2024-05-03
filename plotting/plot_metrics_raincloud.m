% Plot evaluation metrics as rain cloud/half-violin plots. Based on Franziska Pellegrini's fp_plot_pac_univar.m 
% (https://github.com/fpellegrini/PAC/blob/master/plotting/fp_plot_pac_univar.m) function.
%
% Inputs:
%   data         - [integer] (1 x n_iter) data vector
%   color        - [integer] (n_methods x 3) RGB value
%   box_on       - [boolean] whether to show the boxplot 
%   bandwidth    - [float] width of the distribution along the x-axis
%   density_type - ['ks', 'rash'] kernel smoothing function. 'ks' calls MATLAB's ksdensity function, 'rash' uses 'rst_RASH' from Cyril Pernet's
%                  robust statistics toolbox, defauult is 'ks'
%   titles       - [cell array] (1 x n_methods) titles
%   DIROUT       - output directory to save images
%
% Optional input:
%   name - [string] additional string to add to file name

function plot_metrics_raincloud(data, colors, box_on, bandwidth, density_type, titles, DIROUT, varargin)
    
    g = finputcheck(varargin, { ...
            'name'    'string'     { }     '';
            });
    if ischar(g), error(g); end

    figure
%     figone(6,12)

    n_metrics = size(titles, 2);
    o = 1;
    for i_metric = 1:n_metrics

        subplot(1, n_metrics, o)
        
        % raincloud plot with linear y-axis
        h = fp_raincloud_plot_c(data, colors(i_metric, :), box_on, bandwidth, density_type);
    
        view([-90 -90]);
        set(gca, 'Xdir', 'reverse');
        set(gca, 'XLim', [0 1]);  
        htit = title(titles{i_metric});
        htit.Position(1) = -0.12;
        set(gca,'ytick',[])
        ylim([-0.75 2])
        box off

        if o==1
            xlabel('FPR')
            set(gca,'Clipping','Off')
            xt = xticks;
            for ix = xt
                hu = line([ix ix],[2 -15]);
                set(hu, 'color',[0.9 0.9 0.9])
                uistack(hu,'bottom')
            end
            hu1 = line([0 0],[2 -15]);
            set(hu1, 'color',[0 0 0])
        else
            set(gca,'xticklabel',{[]})
            set(gca,'XColor','none','YColor','none','TickDir','out')
            set(gca,'Clipping','Off')
            for ix = xt
                hu = line([ix ix],[2.2 -0.75]);
                set(hu, 'color',[0.9 0.9 0.9])
                uistack(hu,'bottom')
            end
            hu = line([0 0],[2.2 -0.75]);
            set(hu, 'color',[0 0 0])
        end

        o=o+1;
    end
    
    % save figure
    save_err = [DIROUT 'Rain_sim_FPR' g.name '.png'];
    exportgraphics(gcf, save_err)
end