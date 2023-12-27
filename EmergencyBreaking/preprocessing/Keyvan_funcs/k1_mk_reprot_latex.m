function [] = k1_mk_reprot_latex(info,saveDir)
%This functions creates a report from Preprocessind details in ".tex" format
%   Detailed explanation goes here


prepTimeInfo = ['Started on:  ',info.prepStartTime,...
    '\\\\Completed on:  ',info.prepStartTime,'\\\\ \n'];

rawDataInformation = ['Subject ID: ', num2str(info.subject.ID), '\\\\ Eye condition: ', info.subject.eyeCondition,...
    '\\\\ Medical condition: ', info.subject.medicalCondition,' \n',...
    '\\\\ Original sampling frequency: ', num2str(info.data.fs0),' \n',...
    '\\\\ Original data length in mins: ', num2str(info.data.lengthPoints/(info.data.fs0*60)),' \n',...
    '\\\\ Number of channels: ', num2str(info.prep.Nchan0),' \n',...
    '\\\\ Reference electrode: ', info.prep.ref,'\\\\ \n',...
    ];

filteringInformation = ['Low pass filter cut-off: ', num2str(info.prep.lpFilter),'\n',...
    '\\\\ High pass filter cut-off: ', num2str(info.prep.hpFilter),'\n',...
    '\\\\ Band-stop filter Range: [', num2str(info.prep.bsFilter(1)),' ',num2str(info.prep.bsFilter(2)),']\n',...
    '\\\\ Butterworth filter order: ', num2str(info.prep.filterOrder),'\n',...
    ];

downsamplingInformation = ['Downsampling ratio: ', num2str(info.prep.dsRation),'\n',...
    '\\\\ Original sampling frequency: ', num2str(info.data.fs0),'\n',...
    '\\\\ New sampling frequency: ', num2str(info.prep.fs),'\n',...
    ];

epochingInformation = ['Number of epochs: ', num2str(info.prep.Nepochs0),' \n',...
    '\\\\ Epoch length in timeponits: ', num2str(info.prep.epochPoints),' \n',...
    '\\\\ Epoch length in seconds: ', num2str(info.prep.epochPoints/info.prep.fs),' \n',...
    ];

rejectedChanLabels = info.prep.chanLabels(info.outlier.chanRejectFinal');
chanEpochRejection = ['Number of rejected channels: ', num2str(length(find(info.outlier.chanRejectFinal'))),' \n',...
    '\\\\ Indices of rejected channels: ', num2str(find(info.outlier.chanRejectFinal')),' \n',...
    '\\\\ Labels of rejected channels: ', sprintf('%s  ', rejectedChanLabels{:}),' \n',...
    '\\\\ Number of rejected epochs: ', num2str(length(find(info.outlier.epochRejectFinal'))),' \n',...
    '\\\\ Indices of rejected epochs: ', num2str(find(info.outlier.epochRejectFinal)),'\\\\ \n',...
    ];

% ICrejection = ['Number of ICs: ', num2str(length(info.rejectedIC)),' \n',...
%     '\\\\ Number of selected ICs: ', num2str(length(find(~info.rejectedIC))),' \n',...
%     '\\\\ Number of rejected ICs: ', num2str(length(find(info.rejectedIC))),' \n',...
%     '\\\\ Indices of rejected ICs: ', num2str(find(info.rejectedIC)),'\\\\ \n',...
%     ];


chanInterpolation = ['Number of Interpolated channels: ', num2str(length(find(info.outlier.chanRejectFinal'))),' \n',...
    '\\\\ Indices of Interpolated channels: ', num2str(find(info.outlier.chanRejectFinal')),' \n',...
    '\\\\ Labels of Interpolated channels: ', sprintf('%s  ', rejectedChanLabels{:}),'\\\\ \n',...
    ];

spectraPlots = [
   '\\subsection*{ Spectra plots:}\n',...
   '\\subsubsection*{Raw data:}\n',...
   ' \\includegraphics[width=14cm]{spectra_raw.jpg}\\\\\n',...
   '\\subsubsection*{After filterding and downsampling}\n',...
    '\\includegraphics[width=14cm]{spectra_prep1.jpg}\\\\\n',...
    '\\subsubsection*{After Removing extra epochs and Interpolation}\n',...
    '\\includegraphics[width=14cm]{spectra_prep3.jpg}\\\\\n',...
    ];

fileID = fopen([saveDir,'report.tex'],'w');

fprintf(fileID, ['\\documentclass[10pt,a4paper,oneside]{report}\n',...
'\\usepackage{helvet}\n',...
'\\usepackage{helvet}\n',...
'\\usepackage{tgbonum}\n',...
'\\renewcommand{\\familydefault}{\\sfdefault}\n',...
'\\usepackage{array,booktabs}\n',...
'\\usepackage{amsmath}\n',...
'\\usepackage{setspace}\n',...
'\\usepackage{rotating}\n',...
'\\usepackage{color}\n',...
'\\usepackage{multirow}\n',...
'\\usepackage{nth}\n',...
'\\usepackage{graphicx}\n',...
'\\usepackage[font=scriptsize]{caption}\n',...
'\\usepackage[font=scriptsize]{subcaption}\n',...
'\\usepackage{natbib}\n',...
'\\usepackage{geometry}\n',...
'\\usepackage{hyperref}\n',...
'\\geometry{\n',...
	'a4paper,',...
	'total={210mm,297mm},\n',...
	'left=17mm,\n',...
	'right=17mm,\n',...
	'top=20mm,\n',...
	'bottom=20mm,\n',...
'}',...
'\\newcommand*\\rot{\\rotatebox{90}}\n',...
'\\newcommand\\norm[1]{\\left\\lVert#1\\right\\rVert}\n',...
'\\begin{document}\n',...
	'\\author{Keyvan Mahjoory}\n',...
    '\\date{\\taday}\n',...
    '\\section*{Preprocessing Details:}\n',...
    prepTimeInfo,...
    '\\subsection*{1. Data information:}\n',...
    rawDataInformation,...
    '\\subsection*{2. Filtering:}\n',...
    filteringInformation,...
    '\\subsection*{3. Down-sampling:}\n',...
    downsamplingInformation,...
    '\\subsection*{4. Epoching:}\n',...
    epochingInformation,...
    '\\clearpage',...
    '\\subsection*{5. Outlier epoch/channel rejection:}\n',...
    chanEpochRejection,...
    '\\includegraphics[width=14cm]{epoch_chan_deviation.jpg}\\\\\n',...
    '\\includegraphics[width=14cm]{epoch_epoch_deviation.jpg}\\\\\n',...
    '\\includegraphics[width=14cm]{epoch_chan_hfnoise.jpg}\\\\\n',...
    '\\includegraphics[width=14cm]{epoch_epoch_hfnoise.jpg}\\\\\n',...
    '\\includegraphics[width=14cm]{outlier_chans_epochs.jpg}\\\\\n',...
    '\\clearpage',...
   '\\subsection*{7. Channel interpolation:}\n',...
   chanInterpolation,...
    '\\clearpage',...
    spectraPlots,...
    '\\end{document}\n',...
]);
fclose(fileID);


disp(['Output file is saved in:', saveDir]);









end

