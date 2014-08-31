function [fgnumsr1, fgnumsmd, fgnumsfa] = nc_SortByGrowth(coefsPath, fgnums)
% Order fiber groups based on the amount of change
%
% [fgnumsr1, fgnumsmd, fgnumsfa] = nc_SortByGrowth(coefsPath, fgnums)
%
% This function will calculate the amount of change in each qMR parameter 
% (R1, MD, FA, etc.) between childhood and adulthood. The function operates
% on saved model fits and returns a variable denoting the fiber group
% ordering
%
% Inputs:
%
% coefsPath - Path to model saved model fits
% fgnums    - The fiber group numbers (in terms of their order in the AFQ
%             structure), that should be analyzed. I frequently exclude
%             some fiber groups such as the cingulum hippocampus.
% Outputs:
%
% fgnumsr1  - Fiber group numbers ordered in terms of R1 growth
% fgnumsmd  - Fiber group numbers ordered in terms of diffusivity (MD) growth
% fgnumsfa  - Fiber group numbers ordered in terms of RFA growth
%
% Copyright Jason D. Yeatman, August 2014. Code released with:
% Yeatman JD, Wandell BA & Mezer AM (2014). Lifespan maturation 
% and degeneration of human brain white matter. Nature Communications

%% Load model fits
if ~exist('coefsPath','var') || isempty(coefsPath)
    cd(nc_Path)
    load data/coefs_10-Mar-2014.mat
else
    load(coefsPath)
end

%% Set the fiber groups to be the ones used in Yeatman et al. 2014
if ~exist('fgnums','var') || isempty(fgnums)
    fgnums = [1:6 11:28];
end

%% Find the age at which each fiber group reaches "maturity" in terms of R1
% To do this we calculate the vertex of the parabola that we fit to the data

for ii = 1:max(fgnums)
    % vertex equation: x = -b/2a
    v(ii) = -(coefs{3}(1,ii).full(2)./(2*coefs{3}(1,ii).full(1)));
    % Now for each of the fiber groups calculate the amount of change from age
    % 8 to vertex. The variable d will be the difference between the values
    % at these two ages
    d(ii) = diff(polyval(coefs{3}(1,ii).full,[8 v(ii)]));
end

% Remove fiber groups that won't be analyzed
d = d(fgnums);

% Order the fiber groups based on the magnitude of change
[r2,r2i]=sort(d);

% Reorder fgnums based on the magnitude of change. In other words rather
% than considering the ATR the first fiber group, now the fiber group with
% the least change will be the first fiber group.
fgnumsr1 = fgnums(r2i);

%% Find the age at which each fiber group reaches "maturity" in terms of MD
for ii = 1:max(fgnums)
    % equation: x = -b/2a
    vmd(ii) = -(coefs{2}(1,ii).full(2)./(2*coefs{2}(1,ii).full(1)));
    % Now for each of the fiber groups calculate the amount of change from age
    % 8 to vertex
    dmd(ii) = diff(polyval(coefs{2}(1,ii).full,[8 vmd(ii)]));
end

% remove fiber groups that won't be analyzed
dmd = dmd(fgnums);
[r2md,r2imd]=sort(dmd,'descend');
% reorder fgnums based on the magnitude of change
fgnumsmd = fgnums(r2imd);

%% Find the age at which each fiber group reaches "maturity" in terms of FA

for ii = 1:max(fgnums)
    % equation: x = -b/2a
    vfa(ii) = -(coefs{1}(1,ii).full(2)./(2*coefs{1}(1,ii).full(1)));
    % Now for each of the fiber groups calculate the amount of change from age
    % 8 to vertex
    dfa(ii) = diff(polyval(coefs{1}(1,ii).full,[8 vfa(ii)]));
end
% remove fiber groups that won't be analyzed
dfa = dfa(fgnums);
[r2fa,r2ifa]=sort(dfa,'descend');
% reorder fgnums based on the magnitude of change
fgnumsfa = fgnums(r2ifa);