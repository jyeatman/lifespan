function nc_Figure5a

%% Figure 5 compare models
% get age
IdCol = 3; % Column with Westin Havens Ids
Ids = afq.subIds; % Subject Ids
ColNums = [3 7 9 11 13 14 15 23];
[data, header] = organizeBehavioralData(xlsPath, IdCol, Ids, ColNums);
age = data(:,strcmp('Age',header));
% ------------------------ %
% Remove these subjects who had crappy data (03/03/2013
age(excSubs) = nan;
% ------------------------ %

% Loop over fiber tracts and conglomerate data into 1 large vector for each
% tissue property
fgs = AFQ_get(afq,'fgnames')
md_all=[];r1_all=[];tv_all=[];age_all=[];
md_m=[];r1_m=[];
for ii = fgnums
    % Get mean values for this tract
    md_ii = nanmean(AFQ_get(afq,fgs{ii},'md'),2);% - mean(nanmean(AFQ_get(afq,fgs{ii},'md'),2));
    r1_ii = nanmean(AFQ_get(afq,fgs{ii},'R1_2DTI'),2);% - mean(nanmean(AFQ_get(afq,fgs{ii},'R1_2DTI'),2));
    tv_ii = nanmean(AFQ_get(afq,fgs{ii},'TV_map_2DTI'),2);
    
    % concatenate values into 1 large vector
    md_all = vertcat(md_all,md_ii);
    r1_all = vertcat(r1_all,r1_ii);
    tv_all = vertcat(tv_all,tv_ii);
    age_all = vertcat(age_all,age);

    % delete these lines ???
    %     md_m = horzcat(md_m,nanmean(AFQ_get(afq,fgs{ii},'md'),2));
    %     r1_m = horzcat(r1_m,nanmean(AFQ_get(afq,fgs{ii},'R1_2DTI'),2));
end
% remove nans and unwated subject
use = ~isnan(r1_all) & ~isnan(md_all) & ~isnan(age_all);
% ages to calculate values
x0 = 8:.1:80;

% Fit local regression to all the data for all the tracts
r1_curve = localregression(age_all(use),r1_all(use),x0,1,[],25);
md_curve = localregression(age_all(use),md_all(use),x0,1,[],25);
tv_curve = localregression(age_all(use),tv_all(use),x0,1,[],25);

% Fit a second order polynomial to the the local regression predictions
r1_p = polyfit(x0,r1_curve,2);
md_p = polyfit(x0,md_curve,2);
tv_p = polyfit(x0,tv_curve,2);

% Fit a poisson model to the the local regression predictions
r1_pois = fitPoissonCurve(x0,r1_curve);
md_pois = fitPoissonCurve(x0,md_curve);
tv_pois = fitPoissonCurve(x0,tv_curve);

% Bootstrap the local regression model to get confidence intervals
r1_boot = bootstrp(300,@(x,y) localregression(x, y, x0 ,1,[], 25), age_all(use),r1_all(use));
r1_ci = prctile(r1_boot,[16 84]);
md_boot = bootstrp(300,@(x,y) localregression(x, y, x0 ,1,[], 25), age_all(use),md_all(use));
md_ci = prctile(md_boot,[16 84]);
tv_boot = bootstrp(300,@(x,y) localregression(x, y, x0 ,1,[], 25), age_all(use),tv_all(use));
tv_ci = prctile(tv_boot,[16 84]);

% plot r1 curve with poisson and polynomial fits
figure;hold
fill([x0 fliplr(x0)],[r1_ci(1,:) fliplr(r1_ci(2,:))],[.7 .4 .4],'edgecolor',[.7 .4 .4])
plot(x0,polyval(r1_p,x0),'-','linewidth',3,'color',[0 0 0])
plot(x0,evalPoissonCurve(r1_pois,x0),'--','linewidth',3,'color',[0 0 0])
axis tight
xlabel('Age');ylabel('R1 (1/seconds');

% plot md curve with poisson and polynomial fits
figure;hold
fill([x0 fliplr(x0)],[md_ci(1,:) fliplr(md_ci(2,:))],[.4 .4 .7],'edgecolor',[.4 .4 .7])
plot(x0,polyval(md_p,x0),'-','linewidth',3,'color',[0 0 0])
plot(x0,evalPoissonCurve(md_pois,x0),'--','linewidth',3,'color',[0 0 0])
axis tight
xlabel('Age');ylabel('Diffusivity (\mum^2/ms)');

% plot tv curve with poisson and polynomial fits
figure;hold
fill([x0 fliplr(x0)],[tv_ci(1,:) fliplr(tv_ci(2,:))],[.4 .4 .7],'edgecolor',[.4 .4 .7])
plot(x0,polyval(tv_p,x0),'-','linewidth',3,'color',[0 0 0])
plot(x0,evalPoissonCurve(tv_pois,x0),'--','linewidth',3,'color',[0 0 0])
axis tight
xlabel('Age');ylabel('MTV');

% both on same axis
figure; hold
plot(x0,zscore(r1_curve),'-','linewidth',5,'color',[.7 .4 .4])
plot(x0,zscore(-md_curve),'-','linewidth',5,'color',[.4 .4 .7])
%plot(x0,zscore(tv_curve),'-','linewidth',5,'color',[.4 .7 .4])
axis tight
xlabel('Age');
ylabel('Z score');

% R2 across tracts
m = [3 1 5];bar_colors = [.7 .4 .4;.4 .4 .7;.4 .7 .4]
R2_bars = [median(R2{3}(fgnums,m));median(R2{2}(fgnums,m));median(R2{4}(fgnums,m))];
R2_e(1,:) = std(R2{3}(fgnums,m)-repmat(mean(R2{3}(fgnums,m),2),1,3))./sqrt(length(fgnums));
R2_e(2,:) = std(R2{2}(fgnums,m)-repmat(mean(R2{2}(fgnums,m),2),1,3))./sqrt(length(fgnums));
R2_e(3,:) = std(R2{4}(fgnums,m)-repmat(mean(R2{4}(fgnums,m),2),1,3))./sqrt(length(fgnums));

R2_ci(:,:,1) = R2_bars - R2_e;
R2_ci(:,:,2) = R2_bars + R2_e;
errorbargraph2(R2_bars,R2_ci,bar_colors,.3, [], [], [],5)
axis([.5 2.5 15 43]);
ylabel('Cross-validated R^2','fontname','times')
set(gca,'fontname','times','xtick',[1 2],'xticklabel',{'R1' 'Diffusivity'})