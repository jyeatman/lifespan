function nc_Figure4(age_y, coefsPath)
% Figure 4 - Plots comparing lifespan R1 and MD curves
%
% nc_Figure4(age_y, coefsPath)
% 
% Each qMRI parameter takes a unique lifespan trajectory due to the fact
% that there are multple, independent biological processes driving changes
% in the white matter. This function will reproduce the plots in Figure 4
% showing that (a) R1 and diffusivity development rates are independent for
% each tract, (b) mature R1 and diffusivity values reflect different
% properties of the tissue and (c) on average the two parameters take
% different lifespan trajectories (parabola vs. Poisson curve)
%
% Inputs
%
% age_y     - The age to use for childhood. The default is 10 years of age
% coefsPath - Path to the model coefficients
%
%
% Copyright Jason D. Yeatman, August 2014. Code released with:
% Yeatman JD, Wandell BA & Mezer AM (2014). Lifespan maturation 
% and degeneration of human brain white matter. Nature Communications.

% Define the age to consider "childhood"
if ~exist('age_y','var') || isempty(age_y)
    age_y = 10;
end

% Load model fits
if ~exist('coefsPath','var') || isempty(coefsPath)
    cd(nc_Path)
    load data/coefs_10-Mar-2014.mat
else
    load(coefsPath)
end

% Fiber group colors
colors =AFQ_colormap('bgr',24);

% Sort fiber groups by growth and calculate the vertex of each parabola.
% The vertex marks the age of maturity when the function peaks;
[fgnumsr1, fgnumsmd, fgnumsfa, vr1, vmd] = nc_SortByGrowth;

%% Plot change in R1 versus change in diffusivity (MD) (Figure 4a)

% Make a difference of local regression function so that we can calculate
% the amount of change between childhood and maturity
lr_diff = @(x,y,x0) abs(localregression(x,y,x0(2),1,[],20) - localregression(x,y,x0(1),1,[],20));

figure; hold('on')
axis square
c = 0;
for ii = fgnumsr1
    c=c+1;
    % Bootsrap estimates of R1 and diffusivity development (amount of
    % developmental change)
    bs_r1 = bootstrp(1000,@(x,y) lr_diff(x,y,[age_y vr1(ii)]), coefs{3}(3,ii).x,coefs{3}(3,ii).y);
    bs_md = bootstrp(1000,@(x,y) lr_diff(x,y,[age_y vmd(ii)]), coefs{2}(3,ii).x,coefs{2}(3,ii).y);
    
    % Calculate confidence intervals from bootsrapped estimates
    x(c,:) = prctile(bs_r1,[16 50 84]);
    y(c,:) = prctile(bs_md,[16 50 84]);
    
    % Plot
    plot(x(c,[1 3]),[y(c,2) y(c,2)],'-','color',colors(c,:),'linewidth',2);
    plot([x(c,2) x(c,2)],y(c,[1 3]),'-','color',colors(c,:),'linewidth',2);
    plot(x(c,2),y(c,2),'s','color',colors(c,:).*0,'markerfacecolor',colors(c,:).*0,'markersize',4)
end

% Format axes
axis tight
xlabel('R1 developmental change','fontname','times')
ylabel('Diffusivity developmental change','fontname','times')
set(gca,'fontname','times','xtick',.02:.02:.1,'ytick',.02:.01:.1)
axis([0.014    0.0935    0.0244    0.06])
print(sprintf('-f%d',gcf),'-depsc2','-r300','R1_MD')


%% Plot mature R1 versus mature diffusivity (MD) (Figure 4b)

figure; hold('on')
axis square
c = 0;
for ii = fgnumsr1
    c=c+1;
    
    % Bootstrap estimates of mature R1 and mature diffusivity. This is not
    % a difference score, but rather the average value at the functions
    % peak (vertex).
    bs_r1 = bootstrp(1000,@(x,y) localregression(x,y,vr1(ii),1,[],20), coefs{3}(3,ii).x,coefs{3}(3,ii).y);
    bs_md = bootstrp(1000,@(x,y) localregression(x,y,vmd(ii),1,[],20), coefs{2}(3,ii).x,coefs{2}(3,ii).y);
    
    % Calculate confidence intervals
    x(c,:) = prctile(bs_r1,[16 50 84]);
    y(c,:) = prctile(bs_md,[16 50 84]);
    
    % Plot
    plot(x(c,[1 3]),[y(c,2) y(c,2)],'-','color',colors(c,:),'linewidth',2);
    plot([x(c,2) x(c,2)],y(c,[1 3]),'-','color',colors(c,:),'linewidth',2);
    plot(x(c,2),y(c,2),'o','color',colors(c,:).*.6,'markerfacecolor',colors(c,:).*.6,'markersize',8)
end

% Format axes
axis tight
xlabel('Mature R1 (1/seconds)','fontname','times')
ylabel('Mature diffusivity','fontname','times')
set(gca,'fontname','times','xtick',.98:.04:1.14,'ytick',.6:.03:.8)
print(sprintf('-f%d',gcf),'-depsc2','-r300','Mature_R1_MD')

%% Compare average R1 and diffusivity lifespan curves (Figure 4c)

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
x0 = 8:.1:74;

% Fit local regression to all the data for all the tracts
r1_curve = localregression(age_r1,r1_all,x0,1,[],25);
md_curve = localregression(age_md,md_all,x0,1,[],25);

% both on same axis
figure; hold
plot(x0,zscore(r1_curve),'-','linewidth',5,'color',[.7 .4 .4])
plot(x0,zscore(-md_curve),'-','linewidth',5,'color',[.4 .4 .7])
%plot(x0,zscore(tv_curve),'-','linewidth',5,'color',[.4 .7 .4])
axis tight
xlabel('Age');
ylabel('Z score');
grid

