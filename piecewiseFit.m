function [paramsBest, params, sqerr] = piecewiseFit(X,y,c, slope2)
% Fits data with a series of lines joined at hinges
%
% [paramsBest, params, sqerr] = piecewiseFit(X,y,c)
%
% Inputs:
%
% X      - Vector of x measurements
% y      - Vector of y measurements
% c      - Seedpoints for cutpoint or hinge parameters. There are many
%          local minima so it is important to seed the search algorithm 
%          with multiple startpoints. If c is a vector then one hinge will
%          be fit. If c is a matrix with 2 columns than 2 hinges will be
%          fit. In this case each row of the matrix is used as a seed for
%          the two hinges.
% slope2 - Whether or not the slope ofsecond line should be fit or fixed to
%          be flat. The default is for it to be flat. Otherwise slope2 
%          should be set to 'fit'
%
% Outputs:
% paramsBest - The parameters that obtained the lowest error
% params     - The parameters associated with each seedpoint
% sqerr      - The squared error associated with each set of parameters
%
% The parameters are returned as follows
% params(1) = slope 1
% params(2) = intercept 1
% params(3) = slope 2
% params(4) = cutpoint (hinge) 1
% params(5) = slope 3
% params(6) = cutpoint (hinge) 2
%
%
% Copyright Jason D. Yeatman, August 2014. Code released with:
% Yeatman JD, Wandell BA & Mezer AM (2014). Lifespan maturation 
% and degeneration of human brain white matter. Nature Communications.

%% Argument checking
if exist('c','var') && ~isempty(c) && size(c,2) > size(c,1)
    c = c';
end 
if ~exist('slope2','var') || isempty(slope2)
    slope2 = 'flat';
end

%% seed params

if ~exist('c','var') || isempty(c)
    % Estimate the initial slope and interept based on the mean value of X
    c0 = mean(X);
    p0 = polyfit(X(X<c0),y(X<c0),1);
    % Initialize starting parameters
    params0 = [p0 c0];
elseif size(c,2) == 1
    % If an intial cutpoint was passed in then generate a matrix of
    % starting parameters. One set for each defined cutpoint
    for ii = 1:length(c)
        if sum(X<c(ii)) <= 3
            %paramsBest = nan(1,4);
            params0(ii,:) = nan(1,4);
            continue
        end
        % Check whether the second slope should be fit or fixed to zero
%         if strcmp(slope2,'flat')
            params0(ii,:) = [polyfit(X(X<c(ii)),y(X<c(ii)),1) 0 c(ii)];
%         elseif strcmp(slope2,'fit')
%             params0(ii,:) = [polyfit(X(X<c(ii)),y(X<c(ii)),1) polyfit(X(X>=c(ii)),y(X>=c(ii)),1) c(ii)];
%         end
    end
elseif size(c,2) == 2
    [cc1, cc2] = meshgrid(c(:,1),c(:,2));
    c = [cc1(:) cc2(:)];
    for ii = 1:size(c,1)
        if sum(X<c(ii,1)) <= 3 || sum(X>c(ii,2)) <= 3
            %paramsBest = nan(1,6);
            params0(ii,:) = nan(1,6);
            continue
        end
        ptmp = polyfit(X(X>c(ii,2)),y(X>c(ii,2)),1);
        params0(ii,:) = [polyfit(X(X<c(ii,1)),y(X<c(ii,1)),1) 0 c(ii,1) ptmp(1) c(ii,2)];
    end
end

%% Optimization parameters
options = optimset('Display','off');

%% Lower and upper bounds on the parameters
% The second slope term will be fixed to be flat
if size(params0,2) == 4
    paramsLB = [-100 -100 -100 -100];
    paramsUB = [100 100 100 100];
elseif size(params0,2) ==6
    paramsLB = [-100 -100 -100 -100 -100 -100];
    paramsUB = [100 100 100 100 100 100];
end

%% fit piecewise model for each set of starting parameters
for ii = 1:size(params0,1)
    if sum(isnan(params0(ii,:))) > 0
        params(ii,:) = nan;
        sqerr(ii) = nan;
    else
        if size(params0,2) == 4
            if strcmp(slope2,'flat')
                % Fit model with 1 cutpoint and a flat line afer the cut point
                [p, sqerr(ii)] = lsqnonlin(@(params) y - piecewiseEval([params(1:2) 0 params(3)],X),params0(ii,[1 2 4]),paramsLB([1 2 4]),paramsUB([1 2 4]),options);
                % Organize parameters
                params(ii,:) = [p(1) p(2) 0 p(3)];
            elseif strcmp(slope2,'fit')
                % Fit model with 1 cutpoint and a flat line afer the cut point
                [params(ii,:), sqerr(ii)] = lsqnonlin(@(params) y - piecewiseEval(params(1:4),X),params0(ii,:),paramsLB,paramsUB,options);
            end
        elseif size(params0,2) == 6
            if strcmp(slope2,'flat')
               % Fit model with 2 cutpoints
                [p, sqerr(ii)] = lsqnonlin(@(params) y - piecewiseEval([params(1:2) 0 params(3:5)], X),params0(ii,[1 2 4 5 6]),paramsLB([1 2 4 5 6]),paramsUB([1 2 4 5 6]),options);
                % Organize parameters
                params(ii,:) = [p(1) p(2) 0 p(3:5)];
            elseif strcmp(slope2,'fit')
                 % Fit model with 2 cutpoints
                [params(ii,:), sqerr(ii)] = lsqnonlin(@(params) y - piecewiseEval([params(1:6)], X),params0(ii,:),paramsLB,paramsUB,options);
            end
        end
    end
end

%% Choose the best fitting model
[m, idx] = min(sqerr);
paramsBest = params(idx,:);



% figure;hold('on');
% plot(X,y,'o');
% X0 = min(X):1:max(X);
% plot(X0,bilinearEval(params(1),params(2),params(3),params(4),X0),'k-')

return
