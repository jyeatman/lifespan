function [R2, sqErr, yhat, coef] = nc_CrossValidateModels(y, x, model,crossvalidate,bootIter, params)
% Calculate R2 using leave one out cross validation for a variety of models
%
% [R2, sqErr, yhat, coef] = nc_CrossValidateModels(y, x, model,crossvalidate,bootIter, params)
%
% This function will fit and evaluate a number of different types of models 
% using cross validation
%
% The following model classes have been implimented:
% 'linear'          - Linear model
% 'quadratic'       - Second order polynomial
% 'piecewise'       - Piecewise linear model consisting of a line with a
%                     slope joined to a flat line by a hinge
% 'piecewise2'      - Piecewise linear model consisting of 2 lines with
%                     independent slopes each connected by a flat line in 
%                     the middle with 2 independent hinges
% 'piecewisenoflat' - Same as piecewise but the second line has an
%                     independent slope (rather than being flat)
% 'piecewise2noflat'- Same as piecewise2 but the middle line also has a
%                     slope
% 'exponent'        - y = a*x^n + c
% 'lowess'          - local regression model. Requires kendrick kays code.
%                     See https://github.com/knk
% 'poisson'         - A poisson curve
%
% Inputs:
%
% y             - vector of y values
% x             - vector of x values
% model         - string denoting model type (e.g., 'linear')
% crossvalidate - estimate R2 using cross validation? (logical)
% bootIter      - number of bootstrap iterations (scaler >=1)
%
% Outputs:
%
% R2    - Cross validated R2 (model accuracy)
% sqErr - squared error
% yhat  - model prediction
% coef  - model coefficients
%
% Copyright Jason D. Yeatman, August 2014. Code released with:
% Yeatman JD, Wandell BA & Mezer AM (2014). Lifespan maturation 
% and degeneration of human brain white matter. Nature Communications


%% Argument checking
if notDefined('crossvalidate')
    crossvalidate = true;
end
if notDefined('bootIter')
    bootstrap = false;
else
    bootstrap = true;
end

% If a lowess model was requested then parse the bandwidth parameter
if ~isempty(strfind(model,'lowess'))
    if length(model)>6
        params = str2num(model(7:end));
        model = 'lowess';
    end
end
% Set up structure for coeficients
coef = struct('full', [], 'xval', [], 'boot', [],'name',[],'x',[],'y',[]);
% save points in this struct
coef.x = x; coef.y=y;
% And other outputs
R2 = []; sqErr = []; yhat = [];

%% Fit model

% First generate cross validation indices
for ii = 1:length(x)
    ind(:,ii) = horzcat(1:(ii-1), (ii+1):length(x))';
    leftout(ii) = ii;
end

switch(model)
    
    case {'linear' 'lin'}
        coef.name = model;
        % Make regressor matrix by concatenating x with a column of ones
        X = horzcat(x, ones(length(x),1));
        
        % Fit the full model without leaving out data
        coef.full = regress(y, X);
        
        % Bootstrap
        if bootstrap == 1
            coef.boot = bootstrp(bootIter, @regress, y, X);
        end
        
        % Cross validate
        if crossvalidate == 1
            for ii = 1:size(ind,2)
                % Fit the model to the data leaving one out
                coef.xval(:,ii) = regress(y(ind(:,ii)), X(ind(:,ii),:));
                % Predict the y value for the one that was left out
                yhat(ii)   = coef.xval(1,ii) .* X(leftout(ii),1) + coef.xval(2,ii);
            end
        end
    case {'quadratic' 'quad'}
        coef.name = model;
        % Make regressor matrix by concatenating x^2, x and a constant
        X = horzcat(x.^2, x, ones(length(x),1));
        
        % Fit the full model without leaving out data
        coef.full = regress(y, X);
        
        % Bootstrap
        if bootstrap == 1
            coef.boot = bootstrp(bootIter, @regress, y, X);
        end
        
        % Cross validate
        if crossvalidate == 1
            for ii = 1:size(ind,2)
                % Fit the model to the data leaving one out
                coef.xval(:,ii) = regress(y(ind(:,ii)), X(ind(:,ii),:));
                % Predict the y value for the one that was left out
                yhat(ii)   = coef.xval(1,ii) .* X(leftout(ii),1) ...
                    + coef.xval(2,ii) .* X(leftout(ii),2) + coef.xval(3,ii);
            end
        end
        
    case {'piecewise'}
        coef.name = model;
        % We will try cutpoints between age 12 and 30
        c = 12:2:30;
        
        % Fit the full model without leaving out data
        coef.full = piecewiseFit(x, y, c);
        
        % Bootstrap
        if bootstrap == 1
            coef.boot = bootstrp(bootIter,@(x,y) piecewiseFit(x,y,c),x,y);
        end
        
        if crossvalidate == 1
            for ii = 1:size(ind,2)
                % Fit the piecewise model on the subset of the data
                coef.xval(:,ii) = piecewiseFit(x(ind(:,ii)), y(ind(:,ii)), c);
                % Evaluate the model at the left out point
                yhat(ii) = piecewiseEval(coef.xval(:,ii)', x(leftout(ii)));
            end
        end
    case {'piecewise2'}
        coef.name = model;
        % We will try cutpoints between age 12 and 30 and 62 and 80
        %c = [14:1:30; 58:1:74]';
        c = [14:2:28; 60:2:74]';
        % Fit the full model without leaving out data
        coef.full = piecewiseFit(x, y, c);
        
        % Bootstrap
        if bootstrap == 1
            coef.boot = bootstrp(bootIter,@(x,y) piecewiseFit(x,y,c),x,y);
        end
        
        if crossvalidate == 1
            for ii = 1:size(ind,2)
                % Fit the piecewise model on the subset of the data
                coef.xval(:,ii) = piecewiseFit(x(ind(:,ii)), y(ind(:,ii)), c);
                % Evaluate the model at the left out point
                yhat(ii) = piecewiseEval(coef.xval(:,ii)', x(leftout(ii)));
            end
        else
            R2 = nan;
        end
    case {'piecewisenoflat'}
        coef.name = model;
        % We will try cutpoints between age 12 and 30
        c = 12:2:30;
        
        % Fit the full model without leaving out data
        coef.full = piecewiseFit(x, y, c,'fit');
        
        % Bootstrap
        if bootstrap == 1
            coef.boot = bootstrp(bootIter,@(x,y) piecewiseFit(x,y,c,'fit'),x,y);
        end
        
        if crossvalidate == 1
            for ii = 1:size(ind,2)
                % Fit the piecewise model on the subset of the data
                coef.xval(:,ii) = piecewiseFit(x(ind(:,ii)), y(ind(:,ii)), c,'fit');
                % Evaluate the model at the left out point
                yhat(ii) = piecewiseEval(coef.xval(:,ii)', x(leftout(ii)));
            end
        end
    case {'piecewise2noflat'}
        coef.name = model;
        % We will try cutpoints between age 12 and 30 and 62 and 80
        %c = [14:1:30; 58:1:74]';
        c = [14:2:28; 60:2:74]';
        % Fit the full model without leaving out data
        coef.full = piecewiseFit(x, y, c,'fit');
        
        % Bootstrap
        if bootstrap == 1
            coef.boot = bootstrp(bootIter,@(x,y) piecewiseFit(x,y,c,'fit'),x,y);
        end
        
        if crossvalidate == 1
            for ii = 1:size(ind,2)
                % Fit the piecewise model on the subset of the data
                coef.xval(:,ii) = piecewiseFit(x(ind(:,ii)), y(ind(:,ii)), c,'fit');
                % Evaluate the model at the left out point
                yhat(ii) = piecewiseEval(coef.xval(:,ii)', x(leftout(ii)));
            end
        else
            R2 = nan;
        end
    case {'exponent' 'exponential' 'exp'}
        coef.name = model;
        % Set options for nonlinear fitting
        options = optimset('Display','off');
        % Write out the function for the exponential
        expfun = @(p,x) p(1).*x.^p(2) + p(3);
        
        % Do a least squares fit of a line to seed the nonlin fit
        L = polyfit(x, y, 1);
        % Fit the full model without leaving out data
        coef.full = lsqcurvefit(expfun, [L(1) 1 L(2)], x, y, [], [], options);
        
        % Bootstrap
        if bootstrap == 1
            coef.boot = bootstrp(bootIter,@(x,y) lsqcurvefit(expfun, [L(1) 1 L(2)], x, y, [], [], options), x, y);
        end
        
        % Cross validate
        if crossvalidate == 1
            for ii = 1:size(ind,2)
                % Do a least squares fit of a line to seed the nonlin fit
                L = polyfit(x(ind(:,ii)), y(ind(:,ii)), 1);
                
                % Fit the exponential seeding with a second order polynomial
                coef.xval(:,ii) = lsqcurvefit(expfun, [L(1) 1 L(2)], x(ind(:,ii)), y(ind(:,ii)), [], [], options);
                % Evaluate the exponential to get the prediction for the
                % leftout point
                yhat(ii) = feval(expfun, coef.xval(:,ii), x(leftout(ii)));
            end
        end
        
    case {'lowess' 'loes' 'loess'}
        coef.name = model;
        % Transpose x and y into row vectors if necesary
        if size(x,1) > size(x,2)
            x = x';
        end
        if size(y,1) > size(y,2)
            y = y';
        end
        % kernal width
        if notDefined('params')
            k = 15;
        else
            k = params;
        end
        
        % For the coeficients we save first the point at which the model is
        % evaluated and second the predicted value at that point.
        coef.full(:,1) = min(x):.1:max(x);
        coef.full(:,2) = localregression(x, y, coef.full(:,1) ,3,[], k);
        
        % Bootstrap
        if bootstrap == 1
            coef.boot = bootstrp(bootIter,@(x,y) localregression(x, y, coef.full(:,1) ,3,[], k), x, y);
        end
        
        % Crossvalidate
        if crossvalidate == 1
            % Predict y from local regression
            for ii = 1:size(ind,2)
                yhat(ii) = localregression(x(ind(:,ii)), y(ind(:,ii)), x(leftout(ii)) ,[],[],k);
            end
        end
        
    case {'MARS' 'mars'}
        coef.name = model;
        % params structure
        params = aresparams;
        % Linear model not cubic
        params.cubic = false; params.c=3; params.maxFinalFuncs = 3;
        coef.full = aresbuild(x, y,params);
        % Crossvalidate
        if crossvalidate == 1
            % Predict y from local regression
            for ii = 1:size(ind,2)
                model = aresbuild(x(ind(:,ii)), y(ind(:,ii)),params);
                yhat(ii) = arespredict(model, x(leftout(ii)));
            end
        end
        if bootstrap == 1
            coef.boot = bootstrp(bootIter,@(x,y) aresbuild(x,y),x,y)
        end
    case{'poisson'}
        coef.name = model;
        coef.full = fitPoissonCurve(x,y);
        
        % Bootstrap
        if bootstrap == 1
            coef.boot = bootstrp(bootIter, @fitPoissonCurve, x, y);
        end
        
        % Cross validate
        if crossvalidate == 1
            for ii = 1:size(ind,2)
                % Fit the model to the data leaving one out
                coef.xval(:,ii) = fitPoissonCurve(x(ind(:,ii)), y(ind(:,ii),:));
                % Predict the y value for the one that was left out
                yhat(ii) = evalPoissonCurve(coef.xval(:,ii),x(leftout(ii)));
            end
        end
        
        
        
    otherwise
        error('%s is not an implimented model class',model)
end

if crossvalidate == 1
    % match the dimensions of y and yhat
    if size(yhat,2) ~= size(y,2)
        yhat = yhat';
    end
    
    % Calculate the squared error for the model prediction
    sqErr = (y - yhat).^2;
    
    % Calculate crossvalidated R2
    R2 = calccod(yhat,y);
else
    R2 = nan;
end


return