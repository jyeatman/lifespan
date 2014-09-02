function nc_Figure5a(bsIter)
% Figure 5a - Symmetry of R1 over the lifespan
%
% nc_Figure5a([bsIter])
%
% This function will reproduce Figure 5a showing that a parabola is a good
% fit of the R1 lifespan data
%
% Inputs:
% 
% bsIter  - Number of bootstrap iterations
%
% Copyright Jason D. Yeatman, August 2014. Code released with:
% Yeatman JD, Wandell BA & Mezer AM (2014). Lifespan maturation 
% and degeneration of human brain white matter. Nature Communications.

if ~exist('bsIter','var') || isempty(bsIter)
    bsIter = 500;
end
%% Figure 5 compare models

% Get sorted fiber group numbers
[fgnumsr1, fgnumsmd, fgnumsfa] = nc_SortByGrowth;

% Load model fits
if ~exist('coefsPath','var') || isempty(coefsPath)
    cd(nc_Path)
    load data/coefs_10-Mar-2014.mat
else
    load(coefsPath)
end

% Loop over fiber tracts and conglomerate data into 1 large vector for each
% tissue property
md_all=[]; r1_all=[]; age_md=[]; age_r1=[];
for ii = fgnumsr1
    
    % Get mean values for this tract
    md_ii = coefs{2}(1,ii).y;
    r1_ii = coefs{3}(1,ii).y;
    
    % Grab age for each datapoint. There is really no reason we have to
    % grab it from the coefs struct each time, but it seemed easier since a
    % few subjects are missing measurements for a tract
    age_md   = vertcat(age_md,coefs{2}(1,ii).x);
    age_r1   = vertcat(age_r1,coefs{3}(1,ii).x);
    
    % concatenate values into 1 large vector
    md_all = vertcat(md_all, md_ii);
    r1_all = vertcat(r1_all, r1_ii);
end

% ages to calculate values
x0 = 8:.1:80;

% Fit local regression to all the data for all the tracts
r1_curve = localregression(age_r1,r1_all,x0,1,[],25);
md_curve = localregression(age_md,md_all,x0,1,[],25);

% Fit a second order polynomial to the the local regression predictions
r1_p = polyfit(x0,r1_curve,2);
md_p = polyfit(x0,md_curve,2);

% Fit a poisson model to the the local regression predictions
r1_pois = fitPoissonCurve(x0,r1_curve);
md_pois = fitPoissonCurve(x0,md_curve);

% Bootstrap the local regression model to get confidence intervals
r1_boot = bootstrp(bsIter,@(x,y) localregression(x, y, x0 ,1,[], 25), age_r1,r1_all);
r1_ci = prctile(r1_boot,[16 84]);
md_boot = bootstrp(bsIter,@(x,y) localregression(x, y, x0 ,1,[], 25), age_md,md_all);
md_ci = prctile(md_boot,[16 84]);

% Plot r1 curve with poisson and polynomial fits
figure;hold
fill([x0 fliplr(x0)],[r1_ci(1,:) fliplr(r1_ci(2,:))],[.7 .4 .4],'edgecolor',[.7 .4 .4])
plot(x0,polyval(r1_p,x0),'-','linewidth',3,'color',[0 0 0])
plot(x0,evalPoissonCurve(r1_pois,x0),'-','linewidth',3,'color',[.5 .5 .5])
axis tight
xlabel('Age');ylabel('R1 (1/seconds');

% Plot diffusivity curve with poisson and polynomial fits
figure;hold
fill([x0 fliplr(x0)],[md_ci(1,:) fliplr(md_ci(2,:))],[.4 .4 .7],'edgecolor',[.4 .4 .7])
plot(x0,polyval(md_p,x0),'-','linewidth',3,'color',[0 0 0])
plot(x0,evalPoissonCurve(md_pois,x0),'-','linewidth',3,'color',[.5 .5 .5])
axis tight
xlabel('Age');ylabel('Diffusivity (\mum^2/ms)');
