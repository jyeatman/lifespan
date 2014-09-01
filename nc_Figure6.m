function nc_Figure6(coefsPath)
% Code for Figure 6 testing last in first out hypothesis
%
% nc_Figure6(coefsPath)
%
% The cutpoint or hinge parameters (from the piecewise linear model)
% are saved within the coefficients structure. We can loop over the fiber
% groups and plot the age of cutpoint 1 versus cutpoint 2 for each
% fascicle. The last in first out hypothesis predicts that the age of
% maturity (cutpoint 1) will be negatively correlated with the age of
% degeneration (cutpoint 2)
%
% Copyright Jason D. Yeatman, August 2014. Code released with:
% Yeatman JD, Wandell BA & Mezer AM (2014). Lifespan maturation 
% and degeneration of human brain white matter. Nature Communications

% Get the fiber groups in the correct order
[fgnumsr1, fgnumsmd, fgnumsfa] = nc_SortByGrowth;

% Set the colors of each fiber group
colors =AFQ_colormap('bgr',24);

% Load model fits
if ~exist('coefsPath','var') || isempty(coefsPath)
    cd(nc_Path)
    load data/coefs_10-Mar-2014.mat
else
    load(coefsPath)
end

%% Figure 6b -- R1
cut = [];
figure;hold('on');
c=0
for jj = fgnumsr1
    c = c+1;
    if coefs{3}(3,jj).full(4) < coefs{3}(3,jj).full(6)
        x = coefs{3}(2,jj).full(4);
        y = coefs{3}(2,jj).full(6);
        plot(x,y,'ks','markerfacecolor',colors(c,:),'markersize',8)
    end
end
xlabel('Age of maturity');ylabel('Beginning of aging')
title('R1')
axis tight
axis square
set(gcf,'units','normalized','position',[.6 .4 .15 .25],'paperpositionmode','auto')
print(sprintf('-f%d',gcf),'-depsc2','-r300','Figure6b_LIFO.eps')

%% Figure 6c -- Mean diffusivity
cut = [];
figure;hold('on');
c=0
for jj = fgnumsr1
    c = c+1;
    if coefs{2}(3,jj).full(4) < coefs{2}(3,jj).full(6)
        x = coefs{2}(2,jj).full(4);
        y = coefs{2}(2,jj).full(6);
        plot(x,y,'k^','markerfacecolor',colors(c,:),'markersize',8)
    end
end

xlabel('Age of maturity');ylabel('Beginning of aging')
axis tight
axis square
title('Diffusivity')
set(gcf,'units','normalized','position',[.6 .4 .15 .25],'paperpositionmode','auto')
print(sprintf('-f%d',gcf),'-depsc2','-r300','Figure6c_LIFO.eps')


%% Compute correlations between cutpoints
% The last in first out hypothesis predicts that the cutpoints should be
% negatively correlated (i.e., early matureation means later
% degenerations).
cutp_r1 = []; cutp_md = []; cutp_fa = [];
for jj = fgnumsr1
    if R2{3}(jj,2)
        cutp_r1 = vertcat(cutp_r1,[coefs{3}(2,jj).full(4) coefs{3}(2,jj).full(6)]);
    end
    if R2{2}(jj,2)
        cutp_md = vertcat(cutp_md,[coefs{2}(2,jj).full(4) coefs{2}(2,jj).full(6)]);
    end
    cutp_fa = vertcat(cutp_fa,[coefs{1}(2,jj).full(4) coefs{1}(2,jj).full(6)]);
end
disp('R1')
[c p]=corr(cutp_r1)
disp('MD')
[c p]=corr(cutp_md)
disp('FA')
[c p]=corr(cutp_fa)

return