function [R2, coefs, vals_m, age] = nc_ModelSelection(AFQpath, sub_ages, excludeSubs, valName, nodes, showFigs, agerange, bootIter, crossval, outfile,models)
% Test different lifespan models
%
% [R2, coefs, vals_m, age] = nc_ModelSelection(AFQpath, sub_ages, ...
%       excludeSubs, valName, nodes, showFigs, agerange, bootIter, ...
%       crossval, outfile,models)
%
% This code will test a variety of different lifespan models on qMR
% parameters saved in an AFQ structure. Model coefficients will be returned 
% in a strucute that is used by subsequent functions to make figures.
% Currently it is a bit unwealdy and essoteric to the specific analyses of 
% the Nature Communications paper. But it is a good starting place and I 
% will clean it up as time permits. It is mainly a wrapper function for
% nc_CrossValidateModels.m
%
% Inputs
%
% AFQpath     - Path to afq.mat file. See AFQ_run.m
% sub_ages    - A vector of subject ages or a path to an excel files
%               containing subject ages. Must be in the same order as
%               afq.mat
% excludeSubs - Subjects to exclude from model fitting.
% valName     - Name of the value to fit (e.g., 'R1', 'fa', etc.). Must be
%               written as a string and perfectly match how the value is 
%               written in afq.mat
% nodes       - A vector denoting which fiber tract nodes to analyze.
%               Default is nodes = 1:100
% showFigs    - Logical denoting whether figures are desired
% agerange    - Age range over which to fit the model. 
%               Default is [min(sub_ages) max(sub_ages)]
% bootIter    - Number of bootstrap iterations for computing model
%               reliability (scaler)
% crossval    - Whether or not to use leave one out cross validation to
%               estimate model accuracy (logical)
% outfile     - Name of file to save with coefficient estimates
% models      - Cell array of strings denoting which model classes to test.
%               See nc_FitAndEvaluateModels
%
% Outputs:
%
% R2     - Cross validated estimates of model accuracy
% coefs  - Structure of model coefficients
% vals_m - Mean values for fiber tracts
% age    - Ages for each subject
%
% example:
%
% afqPath='/biac4/wandell/biac2/wandell2/data/WH/analysis/AFQ_Callosum_clip_0_02-Sep-2013.mat';
% sub_ages='/biac4/wandell/biac2/wandell2/data/WH/Spreadsheets/Behavorial_Data/BehavioralData_08.19.13_JN.xls';
% [R2, coefs, vals_m, age]= wmdevo_ModelSelection(afqPath, sub_ages, [71 79], 'R1_2DTI', 1:100, 1, [], 100, 1,['coefs' date])
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

%% Get subject ages
% If a vector of ages was not passed in then oad and organize the .xls behavioral spreadsheet

if notDefined('sub_ages')
    sub_ages =  '/biac4/wandell/biac2/wandell2/data/WH/Spreadsheets/Behavorial_Data/BehavioralData_05.14.13_JN.xls';
end

% If rather than a vector of ages, a path to an excel spreadsheet was
% passed in then we can pull the ages out of that. This is essoteric to how
% the data was organized for Yeatman et al., 2014
if ischar(sub_ages) && exist(sub_ages,'file')
    IdCol = 3; % Column with Westin Havens Ids
    Ids = afq.subIds; % Subject Ids
    ColNums = [3 7 9 11 13 14 15 23];
    [data, header] = organizeBehavioralData(sub_ages, IdCol, Ids, ColNums);
    % Find the column for age in the behavioral data
    age = data(:,strcmp('Age',header));
else
    age = sub_ages;
    if length(age) ~= AFQ_get(afq,'number of subjects')
       error('Age must be supplied for each subject \n'); 
    end
end

%% Fit model of lifespan development

% Get fiber group names
fgNames = AFQ_get(afq,'fgnames');

% Exclude subjects if desired
if exist('excludeSubs','var')
    age(excludeSubs) = nan;
end

% Cut age range if desired
if exist('agerange','var') && ~isempty(agerange)
    age(age<agerange(1)) = nan;
    age(age>agerange(2)) = nan;
end

% Loop over fiber groups
for jj = 1:AFQ_get(afq,'numfg')
    fprintf('\nFiting %s:',fgNames{jj});
    % Extract data for the fiber tract
    vals = AFQ_get(afq, fgNames{jj}, valName);
    % Compute the mean for each subject only considering the desired nodes
    % on the tract
    if notDefined('nodewise')
        
        % Calculate the mean value for the desired nodes
        vals_m(:,jj) = nanmean(vals(:,nodes),2);
        
        % Remove nans
        use = ~isnan(age) & ~isnan(vals_m(:,jj));
        X = age(use); y = vals_m(use,jj);
        
        % Try different models and compute error
        for kk = 1:length(models)
            fprintf(', %s',models{kk});
            [R2(jj,kk),~,~,coefs(kk,jj)] = nc_FitAndEvaluateModels(y, X, models{kk},crossval,bootIter);
        end
        % Fit the model for each node instead of the tract mean
    elseif nodewise == 1
        for n = 1:size(vals,2)
            use = ~isnan(age) & ~isnan(vals(:,n));
            X = age(use); y = vals(use,n);
            
            % Try different models and compute error
            for kk = 1:length(models)
                fprintf(', %s',models{kk});
                [R2(jj,kk),~,~,coefs(kk,(jj-1).*100+n)] = nc_FitAndEvaluateModels(y, X, models{kk},crossval,bootIter);
            end
        end
    end
end

% Show a figure if desired
if showFigs == 1
    for kk = 1:length(models)
        wmdevo_PlotModelFits(age, vals_m, coefs(kk,:), models{kk},afq.subIds,valName);
    end
end

% Save output if desired
if exist('outfile','var') && ~isempty(outfile)
    save(outfile);
end
