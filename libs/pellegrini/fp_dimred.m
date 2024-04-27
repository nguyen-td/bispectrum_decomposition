function signal_roi = fp_dimred(signal_sensor,D,A,t)
% Dimensionality reduction within regions of the Desikan-Killiany atlas 
% D is structure resulting from fp_get_Desikan 
% A is inverse filter 
% t is flag that should be set to 1 to select the TRUEVOX pipeline
% Note that here, the number of PCs is fixed to 1. 
%
% Copyright (c) 2023 Franziska Pellegrini and Stefan Haufe

[n_sensors, l_epoch, n_trials] = size(signal_sensor);
ndim = size(A,2);

signal_roi = [];

ipip = 1; %fix number of PCs to 1. 

%loop over regions
for aroi = 1:D.nroi
    
    clear A_ signal_source    
    A_ = A(:, :,D.ind_roi_cortex{aroi},:);
    if t %truevox pipeline
        A_ = A_(:,:,D.sub_ind_roi_region{aroi},:);
    end
    
    %number of voxels at the current roi
    nvoxroi(aroi) = size(A_,3);
    A2{aroi} = reshape(A_, [n_sensors, ndim*nvoxroi(aroi)]);
    
    %project sensor signal to voxels at the current roi (aroi)
    signal_source = A2{aroi}' * signal_sensor(:,:);
    
    %do PCA
    clear signal_roi_ S
    [signal_roi_,S,~] = svd(double(signal_source(:,:))','econ');
%     
%     figure; 
%     pwelch(signal_roi_(:,:),200,[],[],200)
%     legend('1')
%     
    %fixed number of pcs
    npcs(aroi) = ipip;
    
    %bring signal_roi to the shape of npcs x l_epoch x n_trials
    signal_roi = cat(1,signal_roi,reshape((signal_roi_(:,1:npcs(aroi)))',[],l_epoch,n_trials));
    
end