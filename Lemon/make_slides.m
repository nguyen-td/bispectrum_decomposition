function make_slides(png_files, savepath, ppt)
    import mlreportgen.ppt.*

    for p = 1:length(png_files)
        filename = png_files(p).name;
        newslide = add(ppt, 'Blank');
        fig = Picture([savepath, filename]);
        %fig.Width = 'px'; fig.Height = 'px';
        fig.Width = strcat(string(str2double(fig.Width(1:end-2)) * 0.8), 'px'); fig.Height = strcat(string(str2double(fig.Height(1:end-2)) * 0.8), 'px');
        add(newslide, fig);
    end
end