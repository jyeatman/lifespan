function nc_Figure2
% Code for Figure 2 in Yeatman et al. 2014
%
% nc_Figure2
%
% This code will reproduce the montage of plots in Figure 2
%
% Copyright Jason D. Yeatman, August 2014. Code released with:
% Yeatman JD, Wandell BA & Mezer AM (2014). Lifespan maturation 
% and degeneration of human brain white matter. Nature Communications

% Load model fits
if ~exist('coefsPath','var') || isempty(coefsPath)
    cd(nc_Path)
    load data/coefs_10-Mar-2014.mat
else
    load(coefsPath)
end

% Load example fibers and fiber group names
load exampleFibers.mat

% Set the colors and fiber group order
colors =AFQ_colormap('bgr',24);
[fgnumsr1, fgnumsmd, fgnumsfa] = nc_SortByGrowth;

nc_PlotModelFits(coefs{3}(1,:),'R1',fgNames,fgnumsr1,colors);
set(gcf,'units','inches','position',[25 2 14 7],'paperpositionmode','auto');
print(sprintf('-f%d',gcf),'-djpeg','-r300','Figure2_R1.png');
