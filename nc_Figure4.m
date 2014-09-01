function nc_Figure4(age_y, coefsPath)
% Figure 4 - Plots comparing R1 and MD - Multpe process hypothesis

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

