function nc_Figure3
% Code for Figure 3 in Yeatman et al. 2014
%
% nc_Figure3
%
% This code will reproduce plot of predicted versus measured MTV change
% in Figure 3.
%
% Copyright Jason D. Yeatman, August 2014. Code released with:
% Yeatman JD, Wandell BA & Mezer AM (2014). Lifespan maturation 
% and degeneration of human brain white matter. Nature Communications

% Define colors for the points on the plot
colors = AFQ_colormap('bgr',24);

% Sort fiber groups by amount of growth
[fgnumsr1, fgnumsmd, fgnumsfa] = nc_SortByGrowth;

% Load coefficients
cd(nc_Path)
load data/coefs_10-Mar-2014.mat

% For each fiber group calculate the amount of R1 and MTV development
% between childhood and maturity
for ii = 1:length((coefs{3}(1,:)))
    % "Maturity" is at peak R1 definde by the equation for the parabola
    % vertex: x = -b/2a
    v(ii) = -(coefs{3}(1,ii).full(2)./(2*coefs{3}(1,ii).full(1)));
    
    % We use a local regression to estimate the average R1 and MTV values 
    % in childhood (age 10) and maturity. We bootsrap this estimate to get
    % confidence intervals
    bs_r1 = bootstrp(100, @(x,y) localregression(x,y,[10 v(ii)],1,[],20), coefs{3}(3,ii).x,coefs{3}(3,ii).y);
    bs_tv = bootstrp(100, @(x,y) localregression(x,y,[10 v(ii)],1,[],20), coefs{4}(3,ii).x,coefs{4}(3,ii).y);

    % Based on the equation in Mezer et al., 2013 we can predict MTV from
    % R1. We do this for each bootsrap replication of the R1 values
    bs_tvp = 1-1./(bs_r1 .* 0.42 + 0.95);
    
    % Then the difference between the values at maturity and childhood are
    % calculated. Note that these differences are on the bootstrapped
    % replications so we obtain confidence intervals on the difference
    % estimates
    tv_p(ii,:) = prctile(bs_tvp(:,2)-bs_tvp(:,1),[16 50 84]);
    y(ii,:) = prctile(bs_tv(:,2) - bs_tv(:,1),[16 50 84]);   
end

%% Make plots of predicted versus measured MTV change with error bars

figure;hold
c = 0; % start a count

% Loop over the fiber groups
for ii = fgnumsr1
    c = c+1;
    plot(tv_p(ii,:),[y(ii,2) y(ii,2) y(ii,2)],'-','color',colors(c,:),'linewidth',2);
    plot([tv_p(ii,2) tv_p(ii,2) tv_p(ii,2)],y(ii,:),'-','color',colors(c,:),'linewidth',2);
    plot(tv_p(ii,2),y(ii,2),'s','color',colors(c,:).*0,'markerfacecolor',colors(c,:).*0,'markersize',4);
end
axis([0 .021 0 .021]);
axis square
set(gca,'xtick',0:.005:.02,'ytick',0:.005:.02)
xlabel('MTV change predicted from R1')
ylabel('MTV change measured');

return
