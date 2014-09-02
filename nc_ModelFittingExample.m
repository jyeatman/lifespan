% Example script for fitting lifespan curves to qMRI data.
% This script shows how the model coefficients (data/coefs_10-Mar-2014.mat)
% are calculated from data saved in an afq.mat structure.
%
%
% Copyright Jason D. Yeatman, August 2014. Code released with:
% Yeatman JD, Wandell BA & Mezer AM (2014). Lifespan maturation 
% and degeneration of human brain white matter. Nature Communications.


% qMR parameters to analyze (written as they are saved in the AFQ
% structure). To analyze qMRI parameters with AFQ see: 
% http://white.stanford.edu/newlm/index.php/AFQ#Analyzing_quantitative_MRI_data_with_AFQ
valNames = {'fa' 'md' 'R1_2DTI' 'TV_map_2DTI'};

afqPath='path to afq.mat file'; % See AFQ_run.m

% Vector containing each subject's age 
sub_ages = [];

% Number of bootstrap replications to assess model reliability
bootIter = 100;

% Models to fit
models = {'quadratic' 'piecewise2'  'lowess20' 'piecewisenoflat' 'poisson'};

% Loop over the different qMR parameters and fit+evaluate each class of
% model
for ii = 1:length(valNames)
    [R2{ii}, coefs{ii}]= wmdevo_ModelSelection(afqPath, xlsPath, excSubs, valNames{ii}, 1:100, 0, [], bootIter, 1,[],models);
end