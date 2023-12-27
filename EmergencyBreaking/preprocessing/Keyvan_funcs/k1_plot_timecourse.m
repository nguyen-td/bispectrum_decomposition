
%% Developing a (semi)automatic pipeline for epoch/channel rejection

% Implementation Hints
% Adding push buttons to a plot without using GUI
% https://nl.mathworks.com/help/matlab/ref/uicontrol.html
% https://nl.mathworks.com/help/matlab/ref/uipanel.html
% https://nl.mathworks.com/help/matlab/ref/uibuttongroup.html

% How?
% The visualizer should go epoch by epoch, and provide recommendations
% about quality of each epoch (based on e.g. PREP pipeline). It should provide
% additional plots (spectral, topo, ...) which would be helpfull for
% rejecteing/accepting an outlier epoch.


% Applications:
% All connectivity measures that are calculated from epochs
% LRTC which is calculated from 3s-55s min epochs
% sleep stage scoring (requires different ciriterion)!