function [R2, coefs, vals_m, age]= nc_ModelSelection(AFQpath, xlsPath, excludeSubs, valName, nodes, showFigs, agerange, bootIter, crossval, outfile,models)
% Test different lifespan models
%
% [R2, coefs, vals_m, age]= nc_ModelSelection(AFQpath, xlsPath, ...
%       excludeSubs, valName, nodes, showFigs, agerange, bootIter, ...
%       crossval, outfile,models)
%
% This code will test a variety of different lifespan models on qMR
% parameters saved in an AFQ structure. Model coefficients will be returned 
% in a strucute that is used by subsequent functions to make figures.
% Currently it is a bit unwealdy and essoteric to the specific analyses of 
% the Nature Communications paper. But it is a good starting place and I 
% will clean it up as time permits. It is mainly a wrapper function for
% nc_CrossValidateModels
%
% example:
%
% [R2, coefs, vals_m, age] = wmdevo_ModelSelection([], [], [], 'TV_map', 30:70, 1, [])
%
% afqPath='/biac4/wandell/biac2/wandell2/data/WH/analysis/AFQ_Callosum_clip_0_02-Sep-2013.mat';
% xlsPath='/biac4/wandell/biac2/wandell2/data/WH/Spreadsheets/Behavorial_Data/BehavioralData_08.19.13_JN.xls';
% [R2, coefs, vals_m, age]= wmdevo_ModelSelection(afqPath, xlsPath, [71 79], 'R1_2DTI', 1:100, 1, [], 100, 1,['coefs' date])
%
%
% Copyright Jason D. Yeatman, August 2014. Code released with:
% Yeatman JD, Wandell BA & Mezer AM (2014). Lifespan maturation 
% and degeneration of human brain white matter. Nature Communications.

if notDefined('models')
    % These are the model types we'll try
    %models = {'linear' 'piecewise' 'piecewise2' 'quadratic' 'exp' 'lowess'};
    models = {'quadratic' 'piecewise2'  'lowess20' 'piecewisenoflat' 'poisson'};
end

if notDefined('showFigs')
    showFigs = 1;
end
if notDefined('valName')
    valName = 'TV_map';
end
if notDefined('nodes')
    nodes = 30:70;
elseif ischar(nodes) && strcmp(nodes,'individual')
    nodewise = 1;
end
% if notDefined('outfile')
%     outfile = '/biac4/wandell/users/jyeatman/matlab/wmdevo/CoeficientEstimates';
% end
if ~exist('bootIter','var')
    bootIter = 500;
end
if notDefined('crossval')
    crossval = 1;
end
%% Load the save AFQ data structure
if notDefined('AFQpath')
    load /biac4/wandell/biac2/wandell2/data/WH/analysis/AFQ_clip_0_04-Apr-2013.mat;
else
    load(AFQpath);
end

%% Load and organize the .xls behavioral spreadsheet
if notDefined('xlsPath')
    xlsPath =  '/biac4/wandell/biac2/wandell2/data/WH/Spreadsheets/Behavorial_Data/BehavioralData_05.14.13_JN.xls';
end
IdCol = 3; % Column with Westin Havens Ids
Ids = afq.subIds; % Subject Ids
ColNums = [3 7 9 11 13 14 15 23];
[data, header] = organizeBehavioralData(xlsPath, IdCol, Ids, ColNums);

%% Fit piecewise linear functions to fiber tract development

% Get fiber group names
fgNames = AFQ_get(afq,'fgnames');
% Find the column for age in the behavioral data
age = data(:,strcmp('Age',header));
% Exclude subjects if desired
if exist('excludeSubs','var')
    age(excludeSubs) = nan;
end
% Cut age range if desired
if exist('agerange','var') && ~isempty(agerange)
    age(age<agerange(1)) = nan;
    age(age>agerange(2)) = nan;
end
% Generate a set of cutpoints to seed the fitting with
c = [min(age)+5 : 2 : 30];

% Loop over fiber groups
for jj = 1:AFQ_get(afq,'numfg')
    fprintf('\nFiting %s:',fgNames{jj});
    % Extract data for the fiber tract
    vals = AFQ_get(afq, fgNames{jj}, valName);
    % Compute the mean for each subject only considering the desired nodes
    % on the tract
    if notDefined('nodewise')
        vals_m(:,jj) = nanmean(vals(:,nodes),2);
        % Remove nans
        use = ~isnan(age) & ~isnan(vals_m(:,jj));
        X = age(use); y = vals_m(use,jj);
        % Try different models and compute error
        for kk = 1:length(models)
            fprintf(', %s',models{kk});
            [R2(jj,kk),~,~,coefs(kk,jj)] = nc_CrossValidateModels(y, X, models{kk},crossval,bootIter);
        end
        % Fit the model for each node instead of the tract mean
    elseif nodewise == 1
        for n = 1:size(vals,2)
            use = ~isnan(age) & ~isnan(vals(:,n));
            X = age(use); y = vals(use,n);
            
            % Try different models and compute error
            for kk = 1:length(models)
                fprintf(', %s',models{kk});
                [R2(jj,kk),~,~,coefs(kk,(jj-1).*100+n)] = nc_CrossValidateModels(y, X, models{kk},crossval,bootIter);
            end
        end
    end
end

if showFigs == 1
    for kk = 1:length(models)
        wmdevo_PlotModelFits(age, vals_m, coefs(kk,:), models{kk},afq.subIds,valName);
    end
end
if exist('outfile','var') && ~isempty(outfile)
    save(outfile);
end
