function nc_Figure7
% Figure 7 - Compare R1 values for a single subject on different scanners
%
% nc_Figure7
%
% This function will reproduce the plots in Figure 7 demonstrating the
% reliability of R1 measurements made on different scanners. The associated
% afq structures contain data from different scanners, head coils and
% amounts of motion
%
% Copyright Jason D. Yeatman, August 2014. Code released with:
% Yeatman JD, Wandell BA & Mezer AM (2014). Lifespan maturation 
% and degeneration of human brain white matter. Nature Communications.


%% Get T1 vals and run correlations
cd(fullfile(nc_Path,'data','R1_Reliability'));

% Load afq structure with T1 values from different scanners
load afq_31-Mar-2014.mat

% Get fiber group names
fgNames = AFQ_get(afq,'fgnames');

% And fiber group numbers to be analyzed
fgnums = [1:6 11:28];

% Allocate a matrix to fill with t1 values from each fiber tract
t1vals = nan(4,100*length(fgnums));
c = 0;

% For each fiber tract pull out the t1 values for data collected at:
% Stanford CNI with a 32 channel head coil (T1_map_cni32)
% Stanford CNI with a 8 channel head coil (T1_map_cni8)
% Stanford Lucas Center with a 8 channel head coil (T1_map_lc8)
% Hebrew University Siemens 3T with a 8 channel head coil (T1_map_siem8)
for ii = fgnums
    c = c+100;
    t1vals(1,[c-99:c]) = AFQ_get(afq,fgNames{ii},'T1_map_cni32');
    t1vals(2,[c-99:c])  = AFQ_get(afq,fgNames{ii},'T1_map_cni8');
    t1vals(3,[c-99:c])  = AFQ_get(afq,fgNames{ii},'T1_map_lc8');
    t1vals(4,[c-99:c])  = AFQ_get(afq,fgNames{ii},'T1_map_siem8');
end

% Convert from T1 to R1
r1vals = 1./t1vals;

%% Make figures

ax = [.75 1.3 .75 1.3]; % axes size

%% 8 vs 32 channel head coil (Figure 7a)
figure;hold;colormap(1-gray(256))
hist2d(r1vals(1,:),r1vals(2,:),ax(1):.005:ax(2),ax(1):.005:ax(2));
plot(ax(1:2),ax(1:2),'-r')
axis square;axis(ax);
xlabel('R1 (1/seconds) 32 channel head coil');ylabel('R1 (1/seconds) 8 channel head coil')
colorbar
caxis([0 4])
calccod(r1vals(1,:),r1vals(2,:))

%% GE vs Siemens (Figure 7b)
figure;hold;colormap(1-gray(256))
hist2d(r1vals(2,:),r1vals(4,:),ax(1):.005:ax(2),ax(1):.005:ax(2));
plot(ax(1:2),ax(1:2),'-r')
axis square;axis(ax);
xlabel('R1 (1/seconds) GE');ylabel('R1 (1/seconds) Siemens')
colorbar
caxis([0 4])
calccod(r1vals(2,:),r1vals(4,:))

return

%% Motion (Figure 7c) -- not finished yet
cd /biac4/wandell/biac2/wandell2/data/WH/analysis/motion
load afq_04-Mar-2013.mat

t1vals = AFQ_get(afq,'vals','T1_map_lsq');

% convert to R1
r1vals = 1./t1vals;

% Make figure
ax = [.75 1.3 .75 1.3];

% No motion vs motion
figure;hold;colormap(1-gray(256))
hist2d(r1vals(2,:),r1vals(3,:),ax(1):.005:ax(2),ax(1):.005:ax(2));
plot(ax(1:2),ax(1:2),'-r')
axis square;axis(ax);
xlabel('R1 no head motion');ylabel('R1 head motion')
colorbar
caxis([0 4])
calccod(r1vals(2,:),r1vals(3,:))

