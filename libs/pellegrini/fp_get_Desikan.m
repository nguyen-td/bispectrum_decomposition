% Get atlas, voxel, and ROI indices using the Desikan-Killiany atlas.
%
% Input:
%   iReg - number of active voxels per region
%
% Output:
%   D - struct with various fields containing information
%       D.nroi               - [integer] number of ROIs
%       D.nvox               - [integer] number of voxels, aggregated over all ROIs
%       D.ind_cortex         - [array] indices of all voxels, aggregated over all ROIs
%       D.ind_roi_cortex     - [cell array] indices containing the positions of ind_roi{iROI} (region voxels) in ind_cortex, separately for each ROI
%       D.sub_ind_cortex     - [array] indices of active voxels per ROI, the number of active voxels is determined by the iReg input parameter
%       D.sub_ind_roi        - [cell array] same as above but as a cell array 
%       D.sub_ind_roi_region - [cell array] indices containing the positions of sub_ind_roi{iROI} (indices of active voxels) in ind_roi 
%                                           (cell array containing indices of voxels per region)
%       D.ctr_ind_roi        - [cell array] indices of central voxels per ROI (one voxel per ROI)
%       D.ctr_ind_roi_region - [cell array] position of the central voxel in each region
%       D.ctr_ind_cortex     - [array] same as ctr_ind_roi but as an array
%       D.roi2vox            - [array] mapping of voxels to ROIs, contains region assignments for each voxel
%       D.leadfield          - [array] leadfield tensor, loaded from bs_results.mat
%       D.normals            - [array] cortex.VertNormals, normalized Vertices?
%
% Copyright (c) 2022 Franziska Pellegrini and Stefan Haufe
% Documentation by Tien Dung Nguyen

function D = fp_get_Desikan(iReg)

%load leadfield and cortex structures 
load('bs_results.mat')

%number of ROIs in the Desikan-Kiliany Atlas
nroi = length(cortex.Atlas(3).Scouts);

% all voxel inds
ind_cortex = []; %all voxel indices 
ind_roi = {}; %all voxels grouped by region 
ind_roi_cortex = {}; %index of roi voxels in ind_cortex, grouped by region 

%active voxel inds (referring to original indeces) 
sub_ind_cortex = []; %randomly selected active voxels of each region
sub_ind_roi = {}; %randomly selected active voxels of each region, grouped by region
sub_ind_roi_region = {}; %index of active voxel within region 

%central voxel inds (referring to ind_cortex/ ind_roi/ ind_roi_cortex)
ctr_ind_cortex = []; %central voxels of each region
ctr_ind_roi = {}; %central voxels of each region, grouped by region
ctr_ind_roi_region = {}; %index of central voxel within region 

for iROI = 1:nroi
    
    %all voxels
    ind_roi{iROI} = cortex.Atlas(3).Scouts(iROI).Vertices;
    ind_cortex = cat(1, ind_cortex, ind_roi{iROI});
    [~, ind_roi_cortex{iROI}, ~] = intersect(ind_cortex, ind_roi{iROI}); 
    
    %active voxels
    sub_ind_roi{iROI} = cortex.Atlas(3).Scouts(iROI).Vertices(...
        randperm(numel(cortex.Atlas(3).Scouts(iROI).Vertices),iReg)); 
    sub_ind_cortex = cat(1,sub_ind_cortex, sub_ind_roi{iROI});
    for ii = 1:iReg
        sub_ind_roi_region{iROI}(ii) = find(ind_roi{iROI}==sub_ind_roi{iROI}(ii));
    end
    
    %central voxels
    clear pos mid_point
    pos = cortex.Vertices(ind_roi{iROI},:);
    mid_point = mean(pos,1);
    [~,ctr_ind_roi_region{iROI}] = min(eucl(mid_point,pos));
    ctr_ind_roi{iROI} = ind_roi{iROI}(ctr_ind_roi_region{iROI});
    ctr_ind_cortex(iROI) = ctr_ind_roi{iROI};
end

%maps roi indices to voxels 
nvox = length(ind_cortex);
roi2vox = zeros(nvox,1); 
for iroi = 1:nroi 
    roi2vox(ind_cortex(ind_roi_cortex{iroi})) = iroi;
end
roi2vox(roi2vox==0)=[];

%output variable
D.nroi = nroi;
D.nvox = nvox;
D.ind_cortex = ind_cortex;
D.ind_roi_cortex = ind_roi_cortex;
D.sub_ind_cortex = sub_ind_cortex;
D.sub_ind_roi = sub_ind_roi;
D.sub_ind_roi_region = sub_ind_roi_region;
D.ctr_ind_roi = ctr_ind_roi;
D.ctr_ind_roi_region = ctr_ind_roi_region; 
D.ctr_ind_cortex = ctr_ind_cortex; 
D.roi2vox = roi2vox;
D.leadfield = leadfield;
D.normals = cortex.VertNormals;