function main_bsfit_allsubjects(nshuf)
    subjects = 1:1:37;
    exclude = [1, 2, 6, 7, 10, 13, 20, 24, 26, 32, 36]; % 26 subjects
    
    for isub = 1:length(subjects)
        tic
        if ismember(isub, exclude)
            continue
        else
            fprintf('Subject %d  .............................................................................................. \n', isub)
            main_bsfit(nshuf, isub)
        end
        toc
    end
end